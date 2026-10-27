import numpy as np, pandas as pd, pickle
from rdkit import RDLogger
from featurize import featurize
import readacross
RDLogger.DisableLog("rdApp.*")

rf = pickle.load(open("rf_model.pkl", "rb"))

def predict(smiles):
    """Return (mean pIC50, spread, IC50 in nM, conservative IC50 in nM,
    read-across alpha, nearest-neighbor info) for one molecule.

    `mean` is read-across-corrected: blended with the nearest training-set
    analog's measured value when one is close enough (see readacross.py).
    The conservative value (mean + 1 sd pIC50, stronger predicted block) is
    the one to use for safety decisions -- held-out CiPA testing shows the
    raw mean systematically understates potency for the most dangerous
    compounds."""
    fp = featurize(smiles)
    if fp is None:
        return None
    fp_row = fp.reshape(1, -1)
    votes = np.array([t.predict(fp_row)[0] for t in rf.estimators_])
    raw_mean, sd = votes.mean(), votes.std()
    mean, alpha, neighbor = readacross.correct(smiles, raw_mean)
    return mean, sd, 10 ** (9 - mean), 10 ** (9 - (mean + sd)), alpha, neighbor

# --- test on the held-out CiPA drugs ---
held = pd.read_csv("cipa_heldout.csv")
rows = []
for _, r in held.iterrows():
    out = predict(r["Smiles"])
    if out is None:
        continue
    mean, sd, ic50, ic50_cons, alpha, neighbor = out
    n_name, n_sim, n_pic50 = neighbor if neighbor else ("", 0.0, 0.0)
    rows.append({"name": r["name"], "true_pIC50": round(r["pIC50"], 2),
                 "pred_pIC50": round(mean, 2), "sd": round(sd, 2),
                 "error": round(mean - r["pIC50"], 2),
                 "pred_IC50_nM": round(ic50, 1),
                 "conservative_IC50_nM": round(ic50_cons, 1),
                 "readacross_alpha": round(alpha, 2),
                 "nearest_analog": n_name, "analog_sim": round(n_sim, 2),
                 "analog_pIC50": round(n_pic50, 2)})

res = pd.DataFrame(rows).sort_values("true_pIC50", ascending=False)
print(res.to_string(index=False))
print(f"\nHeld-out RMSE: {np.sqrt((res.error ** 2).mean()):.3f}")
res.to_csv("cipa_predictions.csv", index=False)