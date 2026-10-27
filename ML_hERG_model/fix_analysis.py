"""
Three analyses the basic metrics don't cover:

  A. Bias in BOTH directions -- E[pred|true] vs E[true|pred].
     The second is the one that matters for deployment.
  B. Variance-matched rescaling -- does un-shrinking help or hurt?
  C. Threshold tuning on a VALIDATION split, then scored once on test.
     (Tuning on test and reporting test is leakage.)

Reads dataset.npz and rf_model.pkl.
"""
import numpy as np
import pickle
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, roc_curve, roc_auc_score, confusion_matrix
from scipy.stats import pearsonr

d = np.load("dataset.npz")
X, y, tr, te = d["X"], d["y"], d["train"], d["test"]
rf = pickle.load(open("rf_model.pkl", "rb"))

pred_te, y_te = rf.predict(X[te]), y[te]
bar = "=" * 66


def bias_table(truth, prediction, by, label, edges):
    print(f"\n  Binned by {label}:")
    print(f"    {'bin':>12}  {'n':>5}  {'mean true':>10}  {'mean pred':>10}  {'bias':>8}")
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (by >= lo) & (by < hi)
        if m.sum() < 10:
            continue
        print(f"    {lo:5.1f}-{hi:5.1f}  {m.sum():5d}  {truth[m].mean():10.3f}  "
              f"{prediction[m].mean():10.3f}  {(prediction[m]-truth[m]).mean():+8.3f}")


# ---------------------------------------------------- A. both directions
print(f"{bar}\nA. SHRINKAGE: WHICH DIRECTION ARE YOU BINNING?\n{bar}")

r = pearsonr(y_te, pred_te).statistic
print(f"  SD of true values       {y_te.std():.3f}")
print(f"  SD of predictions       {pred_te.std():.3f}")
print(f"  Ratio                   {pred_te.std()/y_te.std():.3f}")
print(f"  Pearson r               {r:.3f}   <- theory says ratio ~= r")
print(f"\n  Predictions are compressed by about {pred_te.std()/y_te.std():.2f}x.")
print("  For an MSE-optimal predictor this is CORRECT, not a defect.")

edges = [0, 4.5, 5.0, 5.5, 6.0, 6.5, 20]
bias_table(y_te, pred_te, y_te, "TRUE value (shows shrinkage -- expected)", edges)
bias_table(y_te, pred_te, pred_te, "PREDICTED value (the deployment question)", edges)

# regression of truth on prediction: slope 1 == calibrated in the useful sense
A = np.vstack([pred_te, np.ones_like(pred_te)]).T
slope, intercept = np.linalg.lstsq(A, y_te, rcond=None)[0]
print(f"\n  Regression of true on predicted: y = {slope:.3f}*pred + {intercept:.3f}")
print(f"  Slope of 1.0 means: when the model says X, the truth averages X.")
if abs(slope - 1.0) < 0.15:
    print("  -> Slope is close to 1: calibrated for decisions. The shrinkage")
    print("     seen when binning by TRUE value does not mislead in use.")
else:
    direction = "under" if slope > 1 else "over"
    print(f"  -> Slope is {slope:.2f}: the model {direction}-states potency even")
    print("     conditional on its own output. A linear recalibration would help.")

# ---------------------------------------------------- B. un-shrinking
print(f"\n{bar}\nB. DOES UN-SHRINKING HELP?\n{bar}")
mu = pred_te.mean()
scale = y_te.std() / pred_te.std()
pred_rescaled = mu + (pred_te - mu) * scale

for name, p in [("original", pred_te), ("variance-matched", pred_rescaled)]:
    rmse = mean_squared_error(y_te, p) ** 0.5
    top = y_te >= 6.5
    print(f"  {name:18s} RMSE {rmse:.3f}   SD {p.std():.3f}   "
          f"bias on potent (>=6.5) {(p[top]-y_te[top]).mean():+.3f}")

print(f"\n  Rescaling by {scale:.2f}x restores the spread and cuts bias on potent")
print("  compounds, at the cost of higher overall RMSE. That is the real")
print("  trade-off: accuracy on average vs. not understating the dangerous tail.")
print("  For cardiac SAFETY the rescaled version is arguably the better choice.")

# ---------------------------------------------------- C. honest threshold
print(f"\n{bar}\nC. DECISION THRESHOLD, TUNED WITHOUT LEAKAGE\n{bar}")
THRESH_TRUE = 5.0

rng = np.random.default_rng(0)
perm = rng.permutation(len(tr))
n_val = len(tr) // 5
val_idx, sub_idx = tr[perm[:n_val]], tr[perm[n_val:]]

print(f"  Refitting on {len(sub_idx)} molecules, tuning on {len(val_idx)} held-out...")
rf2 = RandomForestRegressor(n_estimators=500, min_samples_leaf=2,
                            n_jobs=-1, random_state=0).fit(X[sub_idx], y[sub_idx])

pv, yv = rf2.predict(X[val_idx]), (y[val_idx] > THRESH_TRUE).astype(int)
fpr, tpr, thr = roc_curve(yv, pv)
t_youden = thr[np.argmax(tpr - fpr)]

# also: the threshold giving at least 70% specificity
ok = np.where((1 - fpr) >= 0.70)[0]
t_spec70 = thr[ok[np.argmax(tpr[ok])]] if len(ok) else t_youden

print(f"  Youden-optimal threshold      {t_youden:.3f}")
print(f"  Threshold for >=70% specificity {t_spec70:.3f}")

yt = (y_te > THRESH_TRUE).astype(int)
print(f"\n  Scored ONCE on the untouched test set (n={len(yt)}):")
print(f"    {'threshold':>22}  {'sens':>6}  {'spec':>6}  {'bal acc':>8}  {'MCC':>7}  {'missed':>7}")
for name, t in [("conventional (5.000)", 5.0),
                (f"Youden ({t_youden:.3f})", t_youden),
                (f"spec>=70% ({t_spec70:.3f})", t_spec70)]:
    tn, fp, fn, tp = confusion_matrix(yt, (pred_te > t).astype(int)).ravel()
    sens, spec = tp / (tp + fn), tn / (tn + fp)
    den = np.sqrt(float(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    mcc = (tp*tn - fp*fn) / den if den else float("nan")
    print(f"    {name:>22}  {sens:6.3f}  {spec:6.3f}  {(sens+spec)/2:8.3f}  {mcc:7.3f}  {fn:7d}")

print(f"\n  AUC is unchanged at {roc_auc_score(yt, pred_te):.3f} -- thresholds do not")
print("  alter ranking quality, only where you cut it.")
print(f"\n{bar}")
