"""
Full evaluation of the hERG model.

Four blocks:
  1. Regression metrics on the scaffold-held-out test set
  2. Classification metrics at the 10 uM (pIC50 = 5) blocker threshold
  3. Uncertainty calibration -- does the forest spread mean anything?
  4. Read-across correction -- does blending with the nearest training-set
     analog actually help on the FULL test set, not just a few anecdotes?

Run after train.py. Reads dataset.npz and rf_model.pkl.
"""
import numpy as np
import pandas as pd
import pickle
from sklearn.metrics import (mean_squared_error, mean_absolute_error, r2_score,
                             roc_auc_score, confusion_matrix, average_precision_score)
from scipy.stats import spearmanr, pearsonr
import readacross

d = np.load("dataset.npz")
X, y, tr, te = d["X"], d["y"], d["train"], d["test"]
rf = pickle.load(open("rf_model.pkl", "rb"))

Xte, yte = X[te], y[te]

# Per-tree predictions: (n_trees, n_samples)
print("Polling 500 trees...")
votes = np.stack([t.predict(Xte) for t in rf.estimators_])
pred = votes.mean(axis=0)
sd = votes.std(axis=0)
err = pred - yte
abserr = np.abs(err)

bar = "=" * 62

# ------------------------------------------------- 1. regression
print(f"\n{bar}\n1. REGRESSION  (scaffold-split test set, n={len(yte)})\n{bar}")
rmse = mean_squared_error(yte, pred) ** 0.5
mae = mean_absolute_error(yte, pred)
base_rmse = mean_squared_error(yte, np.full_like(yte, y[tr].mean())) ** 0.5

print(f"  RMSE                  {rmse:.3f}   (typical error {10**rmse:.1f}x in IC50)")
print(f"  MAE                   {mae:.3f}   (median |err| {np.median(abserr):.3f})")
print(f"  R2                    {r2_score(yte, pred):.3f}")
print(f"  Pearson r             {pearsonr(yte, pred).statistic:.3f}")
print(f"  Spearman rho          {spearmanr(yte, pred).statistic:.3f}")
print(f"  Baseline RMSE         {base_rmse:.3f}   (predict training mean)")
print(f"  Improvement           {100*(1-rmse/base_rmse):.1f}%")
print(f"  Test set SD           {yte.std():.3f}   (R2 depends on this)")
print(f"  Within 0.5 log units  {100*np.mean(abserr <= 0.5):.1f}%  (3.2x in IC50)")
print(f"  Within 1.0 log units  {100*np.mean(abserr <= 1.0):.1f}%  (10x in IC50)")

# bias across the range -- the tail problem, quantified
print("\n  Bias by potency bin (negative = under-predicting):")
bins = [(0, 4.5), (4.5, 5.0), (5.0, 5.5), (5.5, 6.0), (6.0, 6.5), (6.5, 20)]
for lo, hi in bins:
    m = (yte >= lo) & (yte < hi)
    if m.sum() >= 5:
        lbl = f"{lo:.1f}-{hi:.1f}" if hi < 20 else f"{lo:.1f}+   "
        print(f"    pIC50 {lbl}  n={m.sum():4d}   mean error {err[m].mean():+.3f}")

# ------------------------------------------------- 2. classification
THRESH = 5.0   # pIC50 5 == IC50 10 uM, the conventional blocker cutoff
print(f"\n{bar}\n2. CLASSIFICATION  (blocker if pIC50 > {THRESH:.0f}, i.e. IC50 < 10 uM)\n{bar}")
yt_c = (yte > THRESH).astype(int)
pr_c = (pred > THRESH).astype(int)

if yt_c.sum() in (0, len(yt_c)):
    print("  Test set is single-class -- skipping.")
else:
    tn, fp, fn, tp = confusion_matrix(yt_c, pr_c).ravel()
    sens = tp / (tp + fn)
    spec = tn / (tn + fp)
    prec = tp / (tp + fp) if (tp + fp) else float("nan")
    mcc_den = np.sqrt(float(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    mcc = (tp*tn - fp*fn) / mcc_den if mcc_den else float("nan")

    print(f"  Blockers in test set  {yt_c.sum()} / {len(yt_c)}  ({100*yt_c.mean():.1f}%)")
    print(f"  AUC-ROC               {roc_auc_score(yt_c, pred):.3f}")
    print(f"  Average precision     {average_precision_score(yt_c, pred):.3f}")
    print(f"  Accuracy              {(tp+tn)/len(yt_c):.3f}")
    print(f"  Sensitivity (recall)  {sens:.3f}   (blockers correctly caught)")
    print(f"  Specificity           {spec:.3f}   (non-blockers correctly cleared)")
    print(f"  Precision             {prec:.3f}")
    print(f"  Balanced accuracy     {(sens+spec)/2:.3f}")
    print(f"  MCC                   {mcc:.3f}   (1 = perfect, 0 = chance)")
    print(f"\n  Confusion matrix:")
    print(f"                   predicted non-blocker   predicted blocker")
    print(f"    true non-block        {tn:6d}              {fp:6d}")
    print(f"    true blocker          {fn:6d}              {tp:6d}")
    print(f"\n  {fn} missed blockers -- the safety-critical error type.")

# ------------------------------------------------- 3. uncertainty
print(f"\n{bar}\n3. UNCERTAINTY CALIBRATION  (does tree spread mean anything?)\n{bar}")
print(f"  sd range              {sd.min():.3f} - {sd.max():.3f}  (mean {sd.mean():.3f})")
rho_u = spearmanr(sd, abserr).statistic
print(f"  Spearman(sd, |err|)   {rho_u:.3f}   <- should be clearly positive")

print("\n  Error by confidence quartile:")
q = np.quantile(sd, [0, .25, .5, .75, 1.0])
for i in range(4):
    m = (sd >= q[i]) & (sd <= q[i+1] if i == 3 else sd < q[i+1])
    tag = ["most confident", "              ", "              ", "least confident"][i]
    print(f"    Q{i+1} {tag}  sd {q[i]:.2f}-{q[i+1]:.2f}   RMSE {np.sqrt((err[m]**2).mean()):.3f}   n={m.sum()}")

print("\n  Coverage (if calibrated, should match the nominal figure):")
for k, nominal in [(1.0, 68.3), (1.645, 90.0), (1.96, 95.0)]:
    cov = 100 * np.mean(abserr <= k * sd)
    flag = "ok" if abs(cov - nominal) < 10 else ("OVERCONFIDENT" if cov < nominal else "underconfident")
    print(f"    +/- {k:.2f} sd   {cov:5.1f}%   (nominal {nominal:.1f}%)   {flag}")

print(f"\n{bar}")
print("Note: raw forest spread is usually overconfident -- it measures")
print("disagreement between trees, not total predictive error. If coverage")
print("is low, scale sd by a constant fitted on the test set and report")
print("both the raw and calibrated version.")

# scaling factor that would fix it
scale = np.sqrt(np.mean((err / np.maximum(sd, 1e-9)) ** 2))
print(f"\n  Scaling sd by {scale:.2f} would give nominal coverage.")
cov_s = 100 * np.mean(abserr <= scale * sd)
print(f"  After scaling: +/-1 sd covers {cov_s:.1f}% (target 68.3%)")

# ------------------------------------------------- 4. read-across correction
print(f"\n{bar}\n4. READ-ACROSS CORRECTION (full test set, n={len(yte)})\n{bar}")
meta = pd.read_csv("herg_model_input.csv").iloc[te].reset_index(drop=True)

corrected = np.empty_like(pred)
alphas = np.empty_like(pred)
for i, (smi, p) in enumerate(zip(meta["Smiles"], pred)):
    corrected[i], alphas[i], _ = readacross.correct(smi, p)

err_c = corrected - yte
abserr_c = np.abs(err_c)
touched = alphas > 0

print(f"  Test molecules with a close-enough analog (alpha>0): "
      f"{touched.sum()} / {len(yte)} ({100*touched.mean():.1f}%)")
print(f"\n  {'':18s} {'RMSE':>8} {'MAE':>8} {'Spearman':>9}")
print(f"  {'raw RF':18s} {np.sqrt((err**2).mean()):8.3f} {abserr.mean():8.3f} "
      f"{spearmanr(yte, pred).statistic:9.3f}")
print(f"  {'read-across':18s} {np.sqrt((err_c**2).mean()):8.3f} {abserr_c.mean():8.3f} "
      f"{spearmanr(yte, corrected).statistic:9.3f}")
if touched.sum() > 5:
    print(f"\n  Among the {touched.sum()} molecules actually touched by read-across:")
    print(f"    raw RF RMSE         {np.sqrt((err[touched]**2).mean()):.3f}")
    print(f"    read-across RMSE    {np.sqrt((err_c[touched]**2).mean()):.3f}")

print("\n  Bias by potency bin, raw RF vs. read-across (negative = under-predicting):")
bins = [(0, 4.5), (4.5, 5.0), (5.0, 5.5), (5.5, 6.0), (6.0, 6.5), (6.5, 20)]
for lo, hi in bins:
    m = (yte >= lo) & (yte < hi)
    if m.sum() >= 5:
        lbl = f"{lo:.1f}-{hi:.1f}" if hi < 20 else f"{lo:.1f}+   "
        print(f"    pIC50 {lbl}  n={m.sum():4d}   raw {err[m].mean():+.3f}   "
              f"read-across {err_c[m].mean():+.3f}")
