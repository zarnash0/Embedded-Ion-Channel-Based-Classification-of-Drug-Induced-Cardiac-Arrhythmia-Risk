# ML hERG Model

The machine-learning half of this project: predicts a drug's hERG/INa/ICaL
channel block from its chemical structure, then feeds that prediction into
the same reduced Priebe-Beuckelmann cardiac cell model used in `MATLAB_sim/`
to simulate the resulting change in action-potential duration (APD90) --
the cellular analog of QT prolongation and Torsades de Pointes risk.

## Pipeline

1. **Data** -- `download.py` pulls hERG IC50 bioactivity records from
   ChEMBL (target CHEMBL240). `clean.py` filters and de-duplicates them,
   recovering censored ("< X nM") measurements that disproportionately
   capture the most potent, highest-risk compounds.
2. **Features** -- `features.py` / `featurize.py` turn each molecule into a
   Morgan fingerprint + physicochemical descriptor vector, with 25 CIPA
   reference drugs held out as a true external test set and the rest split
   by Murcko scaffold (not randomly).
3. **Model** -- `train.py` fits a random forest (+ gradient boosting) to
   predict pIC50, with inverse-frequency sample weighting. `readacross.py`
   corrects predictions for compounds with a close structural analog in the
   training set. `metrics.py` / `cv_analysis.py` / `fix_analysis.py` /
   `label_quality.py` evaluate it in depth (regression, classification,
   uncertainty calibration, repeated scaffold CV, label-noise sources).
4. **Biophysical simulation** -- `pb_model.py` is a Python port of the
   reduced 6-variable Priebe-Beuckelmann model with Hill-equation drug
   block. `run_all_drugs.py` simulates every reference drug at several
   multiples of its clinical Cmax and validates the result against the
   official CiPA risk categories.
5. **Embedded export** -- `export_drugs_h.py` emits `drugs.h`, a C header
   consumed by `drug_mode.c` (and the PIC32 firmware in `Protothreads_1_2`)
   so the same Hill-block logic runs live on the hardware cardiac cell.

## Known limitations

Three CIPA reference drugs (dofetilide, cisapride, terfenadine) remain
genuinely hard to predict even after the fixes above -- one is a true
extrapolation gap (no close training analog), one is a confirmed activity
cliff (a near-identical molecule has very different potency, which a 2D
fingerprint cannot resolve). These are disclosed rather than patched over;
see `metrics.py` and `run_all_drugs.py`'s output for the numbers.

## Regenerating excluded files

`rf_model.pkl` (~140MB) and `herg_browser.csv` (~35MB) are excluded from
version control for size. Regenerate them with:

```
python download.py    # -> herg_raw.csv
python clean.py       # needs herg_browser.csv (re-export from ChEMBL) -> herg_clean.csv
python features.py    # -> dataset.npz, herg_model_input.csv, cipa_heldout.csv
python train.py       # -> rf_model.pkl, gb_model.pkl
```
