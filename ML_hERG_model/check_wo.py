import pandas as pd
from fit_hill import fit_file

BASE = r"C:\Users\User\Downloads\CiPA-Model-Validation-2018\CiPA-Model-Validation-2018\Hill_Fitting\data"
wo, skipped = fit_file(f"{BASE}\\WO_hERG_block.csv", "WO")

fda = pd.read_csv("Li2017_IC50.csv")[["drug", "hERG_IC50", "hERG_h"]]
j = wo.merge(fda, on="drug", how="inner")
j["log_diff"] = (j.IC50 / j.hERG_IC50).apply(lambda x: round(__import__("math").log10(x), 2))

pd.set_option("display.float_format", lambda x: f"{x:,.2f}")
print(j[["drug", "channel", "hERG_IC50", "IC50", "log_diff", "hERG_h", "h"]].to_string(index=False))