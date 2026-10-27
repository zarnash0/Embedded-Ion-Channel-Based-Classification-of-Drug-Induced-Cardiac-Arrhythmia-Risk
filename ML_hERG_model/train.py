import numpy as np, pandas as pd, pickle
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.dummy import DummyRegressor
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import spearmanr

d = np.load("dataset.npz")
X, y, tr, te = d["X"], d["y"], d["train"], d["test"]

# --- inverse-frequency weighting: potent (high pIC50) compounds are rare in
# ChEMBL and an MSE-trained forest shrinks its predictions toward the dense
# middle of the distribution, which is exactly where cardiac safety can't
# afford it (see held-out CiPA errors). Upweight rare bins, sqrt-tempered so
# no single bin dominates the split criterion.
hist, edges = np.histogram(y[tr], bins=20)
bin_idx = np.clip(np.digitize(y[tr], edges[1:-1]), 0, len(hist) - 1)
freq = hist[bin_idx]
sample_weight = 1.0 / np.sqrt(freq)
sample_weight *= len(sample_weight) / sample_weight.sum()

# --- down-weight censored (">"-floor) labels recovered by clean.py: the true
# pIC50 for those rows may be higher than the recorded floor, so trust them
# less than an exact measurement without discarding the signal entirely.
meta = pd.read_csv("herg_model_input.csv")
censored_tr = meta["censored"].values[tr]
sample_weight *= np.where(censored_tr, 0.5, 1.0)

def report(name, pred):
    rmse = mean_squared_error(y[te], pred) ** 0.5
    print(f"{name:16s} RMSE={rmse:.3f}  R2={r2_score(y[te], pred):.3f}  "
          f"Spearman={spearmanr(y[te], pred).statistic:.3f}")

dummy = DummyRegressor().fit(X[tr], y[tr])
print(f"{'Mean baseline':16s} RMSE={mean_squared_error(y[te], dummy.predict(X[te])) ** 0.5:.3f}")

print("Training random forest (a few minutes)...")
rf = RandomForestRegressor(n_estimators=500, min_samples_leaf=2,
                           n_jobs=-1, random_state=0)
rf.fit(X[tr], y[tr], sample_weight=sample_weight)
report("Random forest", rf.predict(X[te]))
pickle.dump(rf, open("rf_model.pkl", "wb"))

print("Training gradient boosting...")
gb = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.08,
                                   min_samples_leaf=10, random_state=0)
gb.fit(X[tr], y[tr], sample_weight=sample_weight)
report("Grad boosting", gb.predict(X[te]))
pickle.dump(gb, open("gb_model.pkl", "wb"))

print("Saved rf_model.pkl and gb_model.pkl")