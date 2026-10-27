import pandas as pd, numpy as np

NOBLOCK = 1e12
RISK = {0: "low", 1: "intermediate", 2: "high"}

d = pd.read_csv("newCiPA.csv")
d = d[d.therapeutic > 0]          # drop flecainide (no Cmax/risk data)

out = pd.DataFrame({
    "drug":      d["drug"],
    "cmax":      d["therapeutic"],
    "risk":      d["CiPA"].map(RISK),
    "hERG_IC50": d["hERGIC50"],       "hERG_h": d["hERGh"],
    "INa_IC50":  d["Peak_sodiumIC50"],"INa_h":  d["Peak_sodiumh"],
    "ICaL_IC50": d["CaLIC50"],        "ICaL_h": d["CaLh"],
})

# zeros in this file mean "no data" => no block
for ch in ["hERG", "INa", "ICaL"]:
    out.loc[out[f"{ch}_IC50"] <= 0, f"{ch}_h"] = 1.0
    out.loc[out[f"{ch}_IC50"] <= 0, f"{ch}_IC50"] = NOBLOCK
    out.loc[out[f"{ch}_h"] <= 0, f"{ch}_h"] = 1.0

out["safety_margin"] = (out.hERG_IC50 / out.cmax).round(1)
out = out.sort_values("safety_margin").reset_index(drop=True)

pd.set_option("display.float_format", lambda x: f"{x:,.1f}")
print(out[["drug", "risk", "cmax", "hERG_IC50", "safety_margin"]].to_string(index=False))
print("\nDrugs:", len(out))
out.to_csv("drug_table.csv", index=False)