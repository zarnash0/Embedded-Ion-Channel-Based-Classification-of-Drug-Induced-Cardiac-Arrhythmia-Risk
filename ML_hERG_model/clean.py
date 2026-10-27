import pandas as pd
import numpy as np

print("Loading...")
cols = ["Molecule ChEMBL ID", "Molecule Name", "Smiles", "Standard Type",
        "Standard Relation", "Standard Value", "Standard Units",
        "pChEMBL Value", "Data Validity Comment", "Target ChEMBL ID"]
df = pd.read_csv("herg_browser.csv", sep=";", usecols=cols, low_memory=False)
print("Raw rows:", len(df))

# keep only hERG IC50 measurements
df = df[df["Target ChEMBL ID"] == "CHEMBL240"]
df = df[df["Standard Type"] == "IC50"]
print("After IC50 + target filter:", len(df))

# drop flagged / unusable rows. Keep exact ("=") AND "<"/"<=" (censored)
# relations: a censored row means the assay hit its floor concentration
# without full block, which disproportionately happens for the MOST potent
# compounds -- discarding it entirely (the old behaviour) silently removes
# exactly the ultra-potent tail the model is worst at.
df = df[df["Data Validity Comment"].isna()]
df["rel"] = df["Standard Relation"].astype(str).str.strip().str.strip("'\"")
df = df[df["rel"].isin(["=", "<", "<="])]
df = df[df["Standard Units"] == "nM"]
df = df.dropna(subset=["Smiles", "Standard Value"])
print("After validity filters:", len(df))

# pIC50: use ChEMBL's value when present, else compute it. For "<"/"<="
# rows this is a FLOOR on potency (true pIC50 may be higher) -- still a
# directionally correct, useful label for the ultra-potent region.
df["Standard Value"] = pd.to_numeric(df["Standard Value"], errors="coerce")
df = df[df["Standard Value"] > 0]
computed = 9 - np.log10(df["Standard Value"])
df["pIC50"] = pd.to_numeric(df["pChEMBL Value"], errors="coerce").fillna(computed)
df["censored"] = df["rel"] != "="

# exact measurements win where both exist; censored floors only fill in
# molecules with no exact measurement at all
exact = df[~df["censored"]]
censored_only = df[df["censored"] & ~df["Smiles"].isin(exact["Smiles"])]
print(f"Exact measurements: {len(exact)} rows, "
      f"{exact['Smiles'].nunique()} molecules")
print(f"Censored-only molecules recovered: {censored_only['Smiles'].nunique()}")

use = pd.concat([exact, censored_only], ignore_index=True)

# one row per molecule (median across repeat measurements)
agg = use.groupby("Smiles", as_index=False).agg(
        pIC50=("pIC50", "median"),
        n=("pIC50", "size"),
        name=("Molecule Name", "first"),
        chembl_id=("Molecule ChEMBL ID", "first"),
        censored=("censored", "all"))
print("Unique molecules:", len(agg), f"(of which censored-only: {agg['censored'].sum()})")

agg.to_csv("herg_clean.csv", index=False)
print(agg["pIC50"].describe())