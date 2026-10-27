import pandas as pd, numpy as np
from scipy.optimize import curve_fit

SRC = r"C:\Users\User\Downloads\CiPA-Model-Validation-2018\CiPA-Model-Validation-2018\Hill_Fitting\data\manual_validationdrug_block.csv"
df = pd.read_csv(SRC)
df = df[df.conc > 0]

def hill(conc, ic50, h):
    return 100.0 / (1.0 + (ic50 / conc) ** h)

rows, skipped = [], []
for (drug, ch), g in df.groupby(["drug", "channel"]):
    if len(g) < 3 or g.block.max() < 15:
        skipped.append(f"{drug}/{ch} (n={len(g)}, maxblock={g.block.max():.0f}%)")
        continue
    try:
        p0 = [g.conc.median(), 1.0]
        p, _ = curve_fit(hill, g.conc.values, g.block.values, p0=p0,
                         bounds=([1e-3, 0.1], [1e12, 5.0]), maxfev=50000)
        resid = g.block.values - hill(g.conc.values, *p)
        rows.append({"drug": drug, "channel": ch,
                     "IC50": float(f"{p[0]:.4g}"), "h": round(p[1], 3),
                     "n": len(g), "rmse_pct": round(np.sqrt((resid**2).mean()), 1),
                     "maxblock": round(g.block.max(), 0)})
    except Exception:
        skipped.append(f"{drug}/{ch} (fit failed)")

fit = pd.DataFrame(rows).sort_values(["drug", "channel"])
print(fit.to_string(index=False))
print("\nSkipped (treated as no block):")
for s in skipped: print("  ", s)

fit.to_csv("validation_ic50_long.csv", index=False)
print("\nSaved validation_ic50_long.csv")