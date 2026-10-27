import pandas as pd, numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupShuffleSplit
from featurize import featurize

RDLogger.DisableLog("rdApp.*")          # silence RDKit parse warnings

df = pd.read_csv("herg_clean.csv")
print("Molecules in:", len(df))

# --- drop the CiPA test drugs so they never train the model ---
CIPA = ["DOFETILIDE","VERAPAMIL","MEXILETINE","QUINIDINE","BEPRIDIL","SOTALOL",
        "IBUTILIDE","VANDETANIB","AZIMILIDE","DISOPYRAMIDE","DOMPERIDONE",
        "DROPERIDOL","PIMOZIDE","CISAPRIDE","ONDANSETRON","TERFENADINE",
        "ASTEMIZOLE","CLARITHROMYCIN","CHLORPROMAZINE","CLOZAPINE",
        "RISPERIDONE","METOPROLOL","NIFEDIPINE","NITRENDIPINE","RANOLAZINE",
        "TAMOXIFEN","DILTIAZEM","LOPERAMIDE"]
names = df["name"].astype(str).str.upper()
held = df[names.isin(CIPA)]
df = df[~names.isin(CIPA)]
print(f"Held out {len(held)} CiPA molecules; {len(df)} remain")
held.to_csv("cipa_heldout.csv", index=False)

# --- Morgan fingerprint + hERG-relevant physicochemical descriptors ---
fps, keep = [], []
for i, s in enumerate(df["Smiles"]):
    fp = featurize(s)
    if fp is not None:
        fps.append(fp); keep.append(i)
df = df.iloc[keep].reset_index(drop=True)
X = np.stack(fps)
y = df["pIC50"].values
print("Feature matrix:", X.shape)

# --- scaffold split ---
def scaffold(s):
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(smiles=s)
    except Exception:
        return s

groups = df["Smiles"].apply(scaffold).values
tr, te = next(GroupShuffleSplit(n_splits=1, test_size=0.2,
                                random_state=0).split(X, y, groups))
print(f"Train: {len(tr)}   Test: {len(te)}   Unique scaffolds: {len(set(groups))}")

np.savez_compressed("dataset.npz", X=X, y=y, train=tr, test=te)
df.to_csv("herg_model_input.csv", index=False)
print("Saved dataset.npz")