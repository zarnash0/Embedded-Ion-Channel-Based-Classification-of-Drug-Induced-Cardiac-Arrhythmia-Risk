"""
Where does the label noise come from, and how much is recoverable?

  A. How much do assay types disagree for the SAME molecule?
  B. How noisy is a single measurement vs an averaged one?
  C. How much data did the '=' filter throw away?
  D. Is there a usable patch-clamp-only subset?

Diagnostic only -- changes nothing. Reads herg_browser.csv.
"""
import numpy as np
import pandas as pd

bar = "=" * 70
SEP = ";"   # change to "," if the parse looks wrong

cols = ["Molecule ChEMBL ID", "Smiles", "Standard Type", "Standard Relation",
        "Standard Value", "Standard Units", "pChEMBL Value",
        "Data Validity Comment", "Assay Type", "BAO Label",
        "Assay Description", "Assay Organism", "Document ChEMBL ID",
        "Target ChEMBL ID"]

print("Loading raw export...")
raw = pd.read_csv("herg_browser.csv", sep=SEP, usecols=lambda c: c in cols,
                  low_memory=False)
raw = raw[raw["Target ChEMBL ID"] == "CHEMBL240"]
print(f"  hERG rows: {len(raw)}")

# ------------------------------------------------------ C. censoring first
print(f"\n{bar}\nC. WHAT THE '=' FILTER DISCARDED\n{bar}")
ic = raw[raw["Standard Type"] == "IC50"].copy()
ic["rel"] = ic["Standard Relation"].astype(str).str.strip().str.strip("'\"")
print(f"  IC50 rows total: {len(ic)}")
for r, n in ic["rel"].value_counts().head(8).items():
    print(f"    relation {r!r:<8} {n:6d}  ({100*n/len(ic):.1f}%)")

cens = ic[ic.rel.isin([">", ">=", "<", "<="])]
cens_mols = set(cens["Smiles"].dropna()) - set(ic[ic.rel == "="]["Smiles"].dropna())
print(f"\n  Molecules ONLY ever reported as censored: {len(cens_mols)}")
print("  These are usable as 'weak blocker' labels for a classifier, or in a")
print("  censored (Tobit-style) regression. Currently discarded entirely.")

other = raw[raw["Standard Type"] != "IC50"]["Standard Type"].value_counts().head(6)
print(f"\n  Other endpoints available (not currently used):")
for t, n in other.items():
    print(f"    {t:<22} {n:6d}")
print("  Ki correlates strongly with IC50 for hERG and could extend coverage")
print("  via multi-task learning or a fitted offset.")

# ------------------------------------------------------ A. assay spread
print(f"\n{bar}\nA. DO ASSAY TYPES AGREE?\n{bar}")
ex = ic[(ic.rel == "=") & (ic["Standard Units"] == "nM")].copy()
ex["Standard Value"] = pd.to_numeric(ex["Standard Value"], errors="coerce")
ex = ex[ex["Standard Value"] > 0].dropna(subset=["Smiles"])
ex["pIC50"] = 9 - np.log10(ex["Standard Value"])

lab = "BAO Label" if "BAO Label" in ex.columns else "Assay Type"
print(f"  Grouping by {lab!r}:\n")
print(f"    {'assay':<42} {'n':>6} {'median pIC50':>13}")
top = ex[lab].value_counts().head(8)
for a, n in top.items():
    print(f"    {str(a)[:42]:<42} {n:6d} {ex[ex[lab]==a].pIC50.median():13.2f}")

# same molecule, different assay type
multi = ex.groupby("Smiles")[lab].nunique()
shared = multi[multi > 1].index
print(f"\n  Molecules measured under >1 assay type: {len(shared)}")
if len(shared) > 20:
    g = ex[ex.Smiles.isin(shared)].groupby(["Smiles", lab]).pIC50.median().reset_index()
    spread = g.groupby("Smiles").pIC50.agg(lambda v: v.max() - v.min())
    print(f"    median cross-assay spread  {spread.median():.2f} log units "
          f"({10**spread.median():.1f}x in IC50)")
    print(f"    75th percentile            {spread.quantile(.75):.2f} "
          f"({10**spread.quantile(.75):.1f}x)")
    print("  This is noise you inherit by pooling methods. It sets a floor on")
    print("  achievable RMSE that no model change can cross.")

# ------------------------------------------------------ B. replicate noise
print(f"\n{bar}\nB. HOW NOISY IS ONE MEASUREMENT?\n{bar}")
rep = ex.groupby("Smiles").pIC50.agg(["count", "std", "median"])
multi_r = rep[rep["count"] >= 2]
print(f"  Molecules with >=2 measurements: {len(multi_r)} of {len(rep)} "
      f"({100*len(multi_r)/len(rep):.0f}%)")
if len(multi_r):
    print(f"  Median within-molecule SD      {multi_r['std'].median():.3f} log units")
    print(f"  Mean within-molecule SD        {multi_r['std'].mean():.3f}")
    print(f"\n  A single measurement therefore carries roughly "
          f"{multi_r['std'].median():.2f} log units of noise.")
    print("  Your model's RMSE cannot meaningfully go below this -- it is the")
    print("  irreducible error in the labels themselves.")

print(f"\n  Measurement count distribution:")
for lo, hi, lbl in [(1, 1, "1 measurement"), (2, 3, "2-3"), (4, 9, "4-9"), (10, 10**6, "10+")]:
    n = ((rep["count"] >= lo) & (rep["count"] <= hi)).sum()
    print(f"    {lbl:<16} {n:6d}  ({100*n/len(rep):.1f}%)")
print("  Weighting the loss by measurement count would let the model trust")
print("  well-replicated labels more than one-off ones.")

# ------------------------------------------------------ D. clean subset
print(f"\n{bar}\nD. IS THERE A CLEANER SUBSET?\n{bar}")
desc = ex["Assay Description"].astype(str).str.lower()
patch = desc.str.contains("patch", na=False)
binding = desc.str.contains("binding|displace|radioligand|astemizole|dofetilide bind", na=False)
flux = desc.str.contains("flux|rubidium|thallium|fluorescen", na=False)
print(f"  Mentions 'patch'         {patch.sum():6d} rows")
print(f"  Looks like binding       {binding.sum():6d} rows")
print(f"  Looks like flux/dye      {flux.sum():6d} rows")
print(f"  Unclassified             {(~(patch|binding|flux)).sum():6d} rows")

sub = ex[patch]
n_mol = sub["Smiles"].nunique()
print(f"\n  Patch-clamp-only molecules: {n_mol}")
if n_mol > 1500:
    print("  -> Large enough to train on. Worth testing: a smaller, cleaner")
    print("     training set often beats a larger, noisier one.")
else:
    print("  -> Probably too small to train on alone. Better used as a")
    print("     high-quality TEST set, or as a feature flag.")

print(f"\n{bar}")
print("Suggested order of work, highest expected payoff first:")
print("  1. Add an assay-type feature, or train patch-clamp-only")
print("  2. Recover censored data as classification labels")
print("  3. Weight the loss by measurement count")
print("  4. Only then: descriptors, tuning, ensembling")
print(bar)
