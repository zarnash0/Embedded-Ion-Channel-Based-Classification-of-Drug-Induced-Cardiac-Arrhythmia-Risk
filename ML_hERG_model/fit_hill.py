import pandas as pd, numpy as np
from scipy.optimize import curve_fit

MIN_MAXBLOCK = 25.0    # below this, IC50 is extrapolation
H_BOUNDS = (0.3, 2.5)  # physiologically plausible Hill coefficients

def hill(conc, ic50, h):
    return 100.0 / (1.0 + (ic50 / conc) ** h)

def fit_file(path, label):
    df = pd.read_csv(path)
    df = df[df.conc > 0]
    rows, skipped = [], []
    for (drug, ch), g in df.groupby(["drug", "channel"]):
        mb = g.block.max()
        if len(g) < 3 or mb < MIN_MAXBLOCK:
            skipped.append((drug, ch, len(g), round(mb, 0)))
            continue
        try:
            p, _ = curve_fit(hill, g.conc.values, g.block.values,
                             p0=[g.conc.median(), 1.0],
                             bounds=([1e-3, H_BOUNDS[0]], [1e12, H_BOUNDS[1]]),
                             maxfev=50000)
            resid = g.block.values - hill(g.conc.values, *p)
            rows.append({"drug": drug, "channel": ch,
                         "IC50": float(f"{p[0]:.4g}"), "h": round(p[1], 3),
                         "n": len(g), "rmse_pct": round(np.sqrt((resid**2).mean()), 1),
                         "maxblock": round(mb, 0),
                         "confident": mb >= 50.0,
                         "set": label})
        except Exception:
            skipped.append((drug, ch, len(g), round(mb, 0)))
    return pd.DataFrame(rows), skipped

if __name__ == "__main__":
    BASE = r"C:\Users\User\Downloads\CiPA-Model-Validation-2018\CiPA-Model-Validation-2018\Hill_Fitting\data"
    val, sv = fit_file(f"{BASE}\\manual_validationdrug_block.csv", "validation")
    tr,  st = fit_file(f"{BASE}\\manual_trainingdrug_block.csv",  "training")
    allfit = pd.concat([tr, val], ignore_index=True)
    allfit.to_csv("hill_fits_all.csv", index=False)
    print(f"Fitted {len(allfit)} drug-channel pairs "
          f"({allfit.confident.sum()} confident, "
          f"{(~allfit.confident).sum()} low-confidence)")
    print(f"Skipped as no-block: {len(sv) + len(st)}")