"""Sanity-check drug block against known pharmacology."""
from pb_model import pb_reduced_params, apply_drug, simulate, ap_features

# (name, cmax_nM, hERG_IC50, hERG_h, INa_IC50, INa_h, ICaL_IC50, ICaL_h)
INF = 1e12
DRUGS = [
    ("control",      0.0,    INF,      1.0, INF,      1.0, INF,     1.0),
    ("dofetilide",   2.0,    4.9,     0.93, 380.5,   0.89, 260.3,  1.16),
    ("sotalol",   14690.0, 110600.0, 0.77, 1.14e9,  0.51, 7.06e6, 0.87),
    ("verapamil",   81.0,   288.0,   0.96, INF,      1.0, 201.8,  1.10),
    ("mexiletine", 4129.0, 28880.0,  0.92, INF,      1.0, 38243.6, 1.03),
    ("quinidine",  3237.0,   992.0,  0.82, 12329.0, 1.49, 51592.3, 0.59),
]

p0 = pb_reduced_params("epi")
base = None
print(f"{'drug':12s} {'dose':>6s} {'APD90':>8s} {'dAPD90%':>9s} "
      f"{'dVdtmax':>9s} {'hERGblk':>8s} {'CaLblk':>8s} {'Nablk':>7s}")

for name, cmax, kh, hh, kn, hn, kc, hc in DRUGS:
    for mult in ([0] if name == "control" else [1, 10, 30]):
        conc = mult * cmax
        p = apply_drug(p0, conc,
                       {"hERG": kh, "INa": kn, "ICaL": kc},
                       {"hERG": hh, "INa": hn, "ICaL": hc})
        t, V = simulate(p, t_end=5000.0)
        f = ap_features(t, V, 4050.0, 5000.0)
        if base is None:
            base = f["APD90"]
        d = 100.0 * (f["APD90"] - base) / base
        b = p["_block"]
        print(f"{name:12s} {mult:5d}x {f['APD90']:8.1f} {d:8.1f}% "
              f"{f['dVdt_max']:9.1f} {100*b['hERG']:7.1f}% "
              f"{100*b['ICaL']:7.1f}% {100*b['INa']:6.1f}%")
