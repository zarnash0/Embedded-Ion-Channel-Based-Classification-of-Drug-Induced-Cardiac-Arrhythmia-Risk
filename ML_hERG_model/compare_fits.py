import pandas as pd, numpy as np

mine = pd.read_csv("hill_fits_all.csv")
mine = mine[mine["set"] == "training"]
fda = pd.read_csv("Li2017_IC50.csv")

rows = []
for ch in ["hERG", "INa", "ICaL", "INaL", "Ito", "IKs", "IK1"]:
    if f"{ch}_IC50" not in fda.columns:
        continue
    f = fda[["drug", f"{ch}_IC50", f"{ch}_h"]].dropna()
    m = mine[mine.channel == ch][["drug", "IC50", "h", "confident"]]
    j = f.merge(m, on="drug")
    for _, r in j.iterrows():
        rows.append({"drug": r["drug"], "channel": ch,
                     "fda_IC50": r[f"{ch}_IC50"], "my_IC50": r["IC50"],
                     "log_diff": round(np.log10(r["IC50"] / r[f"{ch}_IC50"]), 2),
                     "fda_h": r[f"{ch}_h"], "my_h": r["h"],
                     "confident": r["confident"]})

c = pd.DataFrame(rows).sort_values("log_diff", key=abs, ascending=False)
pd.set_option("display.float_format", lambda x: f"{x:,.2f}")
print(c.to_string(index=False))

conf = c[c.confident]
print(f"\nAll pairs:        median |log10 ratio| = {c.log_diff.abs().median():.2f}")
print(f"Confident only:   median |log10 ratio| = {conf.log_diff.abs().median():.2f}  (n={len(conf)})")
c.to_csv("fit_comparison.csv", index=False)