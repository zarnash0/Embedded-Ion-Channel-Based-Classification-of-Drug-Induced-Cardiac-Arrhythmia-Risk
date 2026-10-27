"""
Simulate every drug in drug_table.csv at several multiples of Cmax,
using (a) measured IC50 values and (b) ML-predicted hERG IC50.

The ML path propagates UNCERTAINTY through the simulator instead of
applying a fixed +/-1 sd margin in log-IC50 space: a fixed margin landed on
the steepest part of the Hill block-vs-concentration curve for some drugs
and produced wildly wrong APD90 swings (mexiletine's predicted effect blew
up to +73% vs a true -2.4% when a 1 sd IC50 margin happened to land just
below its Cmax). Instead, 3 representative points (5th/50th/95th percentile
of the read-across-corrected predictive distribution) are each simulated,
and dAPD90 is reported as a median with a [p05, p95] band -- a band that
actually reflects how the simulator responds, not an input-space guess.

Outputs results_all.csv and two figures.
"""
import math
import numpy as np
import pandas as pd
import pickle

from pb_model import pb_reduced_params, apply_drug, simulate, ap_features

DOSES = [1, 4, 10, 30]
T_END = 5000.0          # 5 beats at 1 Hz
BEAT = (4050.0, 5000.0)  # measure the last beat
CELL = "epi"
Q = [0.05, 0.50, 0.95]  # uncertainty band reported for the ML path

# ------------------------------------------------------------------ baseline
p0 = pb_reduced_params(CELL)
t, V = simulate(p0, T_END)
base = ap_features(t, V, *BEAT)
APD90_CTRL = base["APD90"]
print(f"Control: APD90 = {APD90_CTRL:.1f} ms, dV/dt_max = {base['dVdt_max']:.0f} mV/ms\n")

# ------------------------------------------------------- ML predictions (hERG)
def load_ml():
    """Return {drug_name_lower: [pIC50_p05, pIC50_p50, pIC50_p95]} from the
    RF model's predictive distribution (per-tree votes), read-across
    corrected (see readacross.py): the whole vote distribution is shifted
    so it's centered on the corrected estimate rather than the raw mean,
    which preserves the forest's own uncertainty shape instead of assuming
    a symmetric margin."""
    try:
        from rdkit import RDLogger
        from featurize import featurize
        import readacross
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        print("RDKit not available - skipping ML predictions")
        return {}
    try:
        rf = pickle.load(open("rf_model.pkl", "rb"))
        held = pd.read_csv("cipa_heldout.csv")
    except FileNotFoundError as e:
        print(f"Skipping ML predictions ({e})")
        return {}

    out = {}
    for _, r in held.iterrows():
        fp = featurize(r["Smiles"])
        if fp is None:
            continue
        fp_row = fp.reshape(1, -1)
        votes = np.array([tr.predict(fp_row)[0] for tr in rf.estimators_])
        raw_mean = votes.mean()
        corrected_mean, alpha, neighbor = readacross.correct(r["Smiles"], raw_mean)
        shifted = votes + (corrected_mean - raw_mean)
        out[str(r["name"]).lower()] = np.quantile(shifted, Q)
    print(f"Loaded ML predictions for {len(out)} drugs\n")
    return out

ML = load_ml()

# ----------------------------------------------------------------- main loop
def run_sim(cmax, mult, ic50, hills):
    p = apply_drug(p0, mult * cmax, ic50, hills)
    tt, V = simulate(p, T_END)
    f = ap_features(tt, V, *BEAT)
    dapd = 100 * (f["APD90"] - APD90_CTRL) / APD90_CTRL if f["captured"] else np.nan
    return dapd, f, p

drugs = pd.read_csv("drug_table.csv")
rows = []

for _, d in drugs.iterrows():
    name = d["drug"]
    cmax = d["cmax"]
    ic50_meas = {"hERG": d["hERG_IC50"], "INa": d["INa_IC50"], "ICaL": d["ICaL_IC50"]}
    hills     = {"hERG": d["hERG_h"],    "INa": d["INa_h"],    "ICaL": d["ICaL_h"]}

    for mult in DOSES:
        dapd, f, p = run_sim(cmax, mult, ic50_meas, hills)
        rows.append({
            "drug": name, "risk": d["risk"], "variant": "measured",
            "dose_x": mult, "conc_nM": round(mult * cmax, 3),
            "APD90": round(f["APD90"], 1) if f["captured"] else np.nan,
            "dAPD90_pct": round(dapd, 1) if pd.notna(dapd) else np.nan,
            "dAPD90_p05": np.nan, "dAPD90_p95": np.nan,
            "dVdt_max": round(f["dVdt_max"], 1),
            "captured": f["captured"], "frac_captured": float(f["captured"]),
            "hERG_block_pct": round(100 * p["_block"]["hERG"], 1),
            "ICaL_block_pct": round(100 * p["_block"]["ICaL"], 1),
            "INa_block_pct": round(100 * p["_block"]["INa"], 1),
        })

        if name.lower() in ML:
            pic50_q = ML[name.lower()]  # ascending: [p05, p50, p95]
            dapd_q, capt_q, p_q = [], [], []
            for pic50 in pic50_q:
                ic50_ml = {**ic50_meas, "hERG": 10 ** (9 - pic50)}
                dapd_i, f_i, p_i = run_sim(cmax, mult, ic50_ml, hills)
                dapd_q.append(dapd_i)
                capt_q.append(f_i["captured"])
                p_q.append(p_i)
            rows.append({
                "drug": name, "risk": d["risk"], "variant": "predicted",
                "dose_x": mult, "conc_nM": round(mult * cmax, 3),
                "APD90": np.nan,
                "dAPD90_pct": round(dapd_q[1], 1) if pd.notna(dapd_q[1]) else np.nan,
                "dAPD90_p05": round(dapd_q[0], 1) if pd.notna(dapd_q[0]) else np.nan,
                "dAPD90_p95": round(dapd_q[2], 1) if pd.notna(dapd_q[2]) else np.nan,
                "dVdt_max": np.nan,
                "captured": capt_q[1], "frac_captured": float(np.mean(capt_q)),
                "hERG_block_pct": round(100 * p_q[1]["_block"]["hERG"], 1),
                "ICaL_block_pct": round(100 * p_q[1]["_block"]["ICaL"], 1),
                "INa_block_pct": round(100 * p_q[1]["_block"]["INa"], 1),
            })
    done = len(set(r["drug"] for r in rows))
    print(f"  [{done:2d}/{len(drugs)}] {name}")

res = pd.DataFrame(rows)
res.to_csv("results_all.csv", index=False)
print(f"\nSaved results_all.csv ({len(res)} simulations)")

# -------------------------------------------------------------- quick report
m = res[(res.variant == "measured") & (res.dose_x == 1)]
print("\n--- dAPD90 at 1x Cmax (measured IC50), by CiPA risk ---")
for risk in ["high", "intermediate", "low"]:
    g = m[m.risk == risk].sort_values("dAPD90_pct", ascending=False)
    print(f"\n{risk.upper()}")
    for _, r in g.iterrows():
        v = f"{r.dAPD90_pct:+6.1f}%" if pd.notna(r.dAPD90_pct) else "   n/c"
        print(f"   {r.drug:16s} {v}")

ok = m.dropna(subset=["dAPD90_pct"])
if len(ok) > 3:
    from scipy.stats import spearmanr
    rank = {"low": 0, "intermediate": 1, "high": 2}
    rho = spearmanr(ok.dAPD90_pct, ok.risk.map(rank)).statistic
    print(f"\nSpearman(dAPD90 at 1x, CiPA risk) = {rho:.3f}  (n={len(ok)})")

# measured vs ML-predicted median, WITH the propagated uncertainty band
cmp = res[res.variant.isin(["measured", "predicted"]) & (res.dose_x == 1)]
wide = cmp.pivot_table(index="drug", columns="variant",
                        values=["dAPD90_pct", "dAPD90_p05", "dAPD90_p95"])
w = wide.dropna(subset=[("dAPD90_pct", "measured"), ("dAPD90_pct", "predicted")])
if len(w):
    err = (w[("dAPD90_pct", "predicted")] - w[("dAPD90_pct", "measured")]).abs()
    covered = ((w[("dAPD90_pct", "measured")] >= w[("dAPD90_p05", "predicted")]) &
               (w[("dAPD90_pct", "measured")] <= w[("dAPD90_p95", "predicted")]))
    print(f"\nMeasured vs ML-predicted (median) dAPD90 at 1x Cmax (n={len(w)}):")
    print(f"   mean |error| = {err.mean():.1f} percentage points")
    print(f"   max  |error| = {err.max():.1f}  ({err.idxmax()})")
    print(f"   measured value falls within predicted [p05,p95] band: "
          f"{covered.sum()}/{len(w)} ({100*covered.mean():.0f}%)")
