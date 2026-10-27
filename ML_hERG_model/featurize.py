"""
Shared featurization for the hERG model: Morgan fingerprint plus a handful of
physicochemical descriptors known to drive hERG liability (lipophilicity,
basic nitrogen, aromatic ring count -- the classic hERG pharmacophore).

Used by features.py (training), predict.py and run_all_drugs.py (inference).
All three must featurize identically or predictions silently misalign with
the columns the model was trained on.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, Descriptors, rdMolDescriptors

_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

DESCRIPTORS = [
    ("MolLogP", Descriptors.MolLogP),
    ("TPSA", Descriptors.TPSA),
    ("MolWt", Descriptors.MolWt),
    ("NumAromaticRings", rdMolDescriptors.CalcNumAromaticRings),
    ("NumHAcceptors", Descriptors.NumHAcceptors),
    ("NumHDonors", Descriptors.NumHDonors),
    ("NumRotatableBonds", Descriptors.NumRotatableBonds),
    ("FractionCSP3", rdMolDescriptors.CalcFractionCSP3),
    ("BasicNCount", lambda m: Descriptors.fr_NH0(m) + Descriptors.fr_NH1(m) + Descriptors.fr_NH2(m)),
    ("AromaticNCount", Descriptors.fr_Ar_N),
]
N_DESC = len(DESCRIPTORS)
FP_SIZE = 2048


def featurize(smiles):
    """Morgan fingerprint (2048 bits) + N_DESC physicochemical descriptors.
    Returns a 1D float32 array, or None if the SMILES doesn't parse."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = np.array(_gen.GetFingerprint(mol), dtype=np.float32)
    desc = np.array([f(mol) for _, f in DESCRIPTORS], dtype=np.float32)
    return np.concatenate([fp, desc])
