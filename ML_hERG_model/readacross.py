"""
Read-across correction: blend an RF point prediction with nearest
structural analogs already in the training set.

Why: the RF's 2048-bit fingerprint signal for a single rare, highly potent
analog gets diluted across 500 trees averaging over mostly dissimilar
neighbors. When a close analog genuinely exists (e.g. dofetilide has a
training-set neighbor at Tanimoto 0.61 with pIC50 8.29, vs. dofetilide's
true 7.85 -- a good match), its *measured* value is more informative than
the RF's ensemble average. When no close analog exists (e.g. ibutilide,
max similarity 0.41), there is nothing to read across from and the RF
prediction is left alone.

alpha interpolates linearly between SIM_LOW (pure RF) and SIM_HIGH (pure
read-across) based on the single nearest neighbor's similarity; the k
nearest neighbors above SIM_LOW are then averaged, weighted by similarity,
to produce the read-across estimate itself.
"""
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs, rdFingerprintGenerator

_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

SIM_LOW, SIM_HIGH = 0.30, 0.70
K_NEIGHBORS = 3

_pool = None  # (fingerprints, pIC50 array, name array), lazily loaded


def _load_pool():
    global _pool
    if _pool is not None:
        return _pool
    meta = pd.read_csv("herg_model_input.csv")
    d = np.load("dataset.npz")
    tr = d["train"]
    meta = meta.iloc[tr].reset_index(drop=True)
    fps = [_gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in meta["Smiles"]]
    labels = meta["name"].fillna(meta["chembl_id"]).values
    _pool = (fps, meta["pIC50"].values, labels)
    return _pool


def correct(smiles, model_pred):
    """Return (corrected_pIC50, alpha, neighbor_info).

    alpha is how much weight went to the read-across estimate (0 = pure
    model, 1 = pure read-across). neighbor_info is (name, similarity,
    pIC50) of the single nearest analog, for reporting.
    """
    fps, y_pool, names = _load_pool()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return model_pred, 0.0, None

    q_fp = _gen.GetFingerprint(mol)
    sims = np.array(DataStructs.BulkTanimotoSimilarity(q_fp, fps))

    order = np.argsort(sims)[::-1]
    best_i = order[0]
    neighbor_info = (names[best_i], float(sims[best_i]), float(y_pool[best_i]))

    cand = order[sims[order] >= SIM_LOW][:K_NEIGHBORS]
    if len(cand) == 0:
        return model_pred, 0.0, neighbor_info

    max_sim = sims[cand[0]]
    alpha = float(np.clip((max_sim - SIM_LOW) / (SIM_HIGH - SIM_LOW), 0.0, 1.0))
    if alpha == 0.0:
        return model_pred, 0.0, neighbor_info

    knn_pred = np.average(y_pool[cand], weights=sims[cand])
    corrected = alpha * knn_pred + (1 - alpha) * model_pred
    return float(corrected), alpha, neighbor_info
