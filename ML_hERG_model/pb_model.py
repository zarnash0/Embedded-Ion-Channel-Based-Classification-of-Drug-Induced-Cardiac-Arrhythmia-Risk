"""
Python port of the Bernus et al. (2002) 6-variable reduced
Priebe-Beuckelmann ventricular cell model, with multi-channel drug block.

Ported from the MATLAB implementation:
  pb_constants.m, pb_reduced_params.m, pb_reduced_gates.m,
  pb_reduced_currents.m, simulate_PB_reduced.m

State: V, m, v, f, to, X
Integration: Rush-Larsen exponential update (same scheme as the MATLAB
"exponential" method), dt = 0.02 ms by default.
"""

import math
import numpy as np

# ---------------------------------------------------------------- constants
R = 8.3143        # J/(mol*K)
T = 310.15        # K (37 C)
F = 96.4867       # C/mmol
RTF = R * T / F   # ~26.73 mV


# ------------------------------------------------------------------- params
def pb_reduced_params(cell_type="epi"):
    p = {
        "cellType": cell_type,
        "gNa": 16.0, "gCa": 0.064, "gNab": 0.001, "gCab": 0.00085,
        "gNaK": 1.3, "gNaCa": 1000.0, "gK1": 3.9,
        "Nai": 10.0, "Nae": 138.0, "Ki": 140.0, "Ke": 4.0,
        "Cai": 0.0004, "Cae": 2.0, "fCa": 0.6,
        "KmNai_NaK": 10.0, "KmKo_NaK": 1.5,
        "KmNa_NaCa": 87.5, "KmCa_NaCa": 1.38,
        "ksat_NaCa": 0.1, "eta_NaCa": 0.35,
        "V0": -90.2,
    }
    if cell_type == "epi":
        p.update(gto=0.4, gK=0.019, invP=1.0, Vshift=0.0)
    elif cell_type == "M":
        p.update(gto=0.35, gK=0.013, invP=1.7, Vshift=-4.0)
    elif cell_type == "endo":
        p.update(gto=0.13, gK=0.019, invP=2.8, Vshift=-12.0)
    else:
        raise ValueError('cell_type must be "epi", "M", or "endo"')
    return p


# -------------------------------------------------------------------- gates
def pb_reduced_gates(V, p):
    xm = V + 47.13
    if abs(xm) < 1e-6:
        am = 3.2
    else:
        am = 0.32 * xm / (1.0 - math.exp(-0.1 * xm))
    bm = 0.08 * math.exp(-V / 11.0)
    minf = am / (am + bm)
    taum = 1.0 / (am + bm)

    vinf = 0.5 * (1.0 - math.tanh(7.74 + 0.12 * V))
    tauv = 0.25 + 2.24 * ((1.0 - math.tanh(7.74 + 0.12 * V)) /
                          (1.0 - math.tanh(0.07 * (V + 92.4))))

    ad = 14.9859 * math.exp(-0.5 * ((V - 22.36) / 16.6813) ** 2) / (16.6813 * math.sqrt(2 * math.pi))
    bd = 0.1471 - (5.3 * math.exp(-0.5 * ((V - 6.2744) / 14.93) ** 2)) / (14.93 * math.sqrt(2 * math.pi))
    dinf = ad / (ad + bd)

    af = 6.872e-3 / (1.0 + math.exp((V - 6.1546) / 6.1223))
    bf = (0.0687 * math.exp(-0.1081 * (V + 9.8255)) + 0.0112) / \
         (1.0 + math.exp(-0.2779 * (V + 9.8255))) + 5.474e-4
    finf = af / (af + bf)
    tauf = 1.0 / (af + bf)

    ar = (0.5266 * math.exp(-0.0166 * (V - 42.2912))) / (1.0 + math.exp(-0.0943 * (V - 42.2912)))
    br = (5.186e-5 * V + 0.5149 * math.exp(-0.1344 * (V - 5.0027))) / \
         (1.0 + math.exp(-0.1348 * (V - 5.186e-5)))
    rinf = ar / (ar + br)

    Vs = V - p["Vshift"]
    ato = (5.612e-5 * Vs + 0.0721 * math.exp(-0.173 * (Vs + 34.2531))) / \
          (1.0 + math.exp(-0.1732 * (Vs + 34.2531)))
    bto = (1.215e-4 * Vs + 0.0767 * math.exp(-1.66e-9 * (Vs + 34.0235))) / \
          (1.0 + math.exp(-0.1604 * (Vs + 34.0235)))
    toinf = ato / (ato + bto)
    tauto = p["invP"] / (ato + bto)

    if p["cellType"] == "M":
        Xinf = 0.972 / (1.0 + math.exp(-2.036 - 0.0834 * V))
        tauX = 380.0 * math.exp(-((25.5 + V) ** 2) / 156.0) + \
               166.0 * (1.0 + math.tanh(0.558 + 0.0169 * V))
    else:
        Xinf = 0.988 / (1.0 + math.exp(-0.861 - 0.0620 * V))
        tauXprime = 40.0 * (1.0 - math.tanh(160.0 + 2.0 * V))
        tauX = 240.0 * math.exp(-((25.5 + V) ** 2) / 156.0) + \
               182.0 * (1.0 + math.tanh(0.154 + 0.0116 * V)) + tauXprime

    return dict(minf=minf, taum=taum, vinf=vinf, tauv=tauv,
                dinf=dinf, finf=finf, tauf=tauf, rinf=rinf,
                toinf=toinf, tauto=tauto, Xinf=Xinf, tauX=tauX)


# ------------------------------------------------------------------ Hill law
def unblocked_fraction(conc, ic50, hill=1.0):
    """Fraction of channel still conducting at drug concentration `conc`."""
    if conc <= 0 or ic50 <= 0:
        return 1.0
    return 1.0 / (1.0 + (conc / ic50) ** hill)


def apply_drug(p, conc_nM, ic50s, hills=None):
    """Return a copy of p with gNa/gCa/gK scaled by Hill block.

    ic50s: dict with keys 'hERG', 'INa', 'ICaL' (nM). Missing -> no block.
    hills: dict with same keys (default 1.0 each).
    """
    hills = hills or {}
    q = dict(p)
    mapping = {"hERG": "gK", "INa": "gNa", "ICaL": "gCa"}
    block = {}
    for ch, gname in mapping.items():
        ic50 = ic50s.get(ch, float("inf"))
        h = hills.get(ch, 1.0)
        frac = unblocked_fraction(conc_nM, ic50, h)
        q[gname] = p[gname] * frac
        block[ch] = 1.0 - frac
    q["_block"] = block
    return q


# ---------------------------------------------------------------- simulation
def simulate(p, t_end, dt=0.02, period=1000.0, stim_dur=3.0,
             stim_amp=15.0, stim_offset=50.0, V0=None):
    """Rush-Larsen integration. Returns (t, V) as numpy arrays."""
    n = int(round(t_end / dt)) + 1
    V_out = np.empty(n)

    V = p["V0"] if V0 is None else V0
    g = pb_reduced_gates(V, p)
    m, v, f, to, X = g["minf"], g["vinf"], g["finf"], g["toinf"], g["Xinf"]

    ENa = RTF * math.log(p["Nae"] / p["Nai"])
    ECa = RTF / 2 * math.log(p["Cae"] / p["Cai"])
    Eto = RTF * math.log((0.043 * p["Nae"] + p["Ke"]) / (0.043 * p["Nai"] + p["Ki"]))
    EK = RTF * math.log(p["Ke"] / p["Ki"])

    sigma = (1.0 / 7.0) * (math.exp(p["Nae"] / 67.3) - 1.0)
    fNaK2 = 1.0 / (1.0 + (p["KmNai_NaK"] / p["Nai"]) ** 1.5) * (p["Ke"] / (p["Ke"] + p["KmKo_NaK"]))
    naca_pref = (p["KmNa_NaCa"] ** 3 + p["Nae"] ** 3) ** -1 * (p["KmCa_NaCa"] + p["Cae"]) ** -1

    gNa, gCa, gto, gK, gK1 = p["gNa"], p["gCa"], p["gto"], p["gK"], p["gK1"]
    gNab, gCab, gNaK, gNaCa, fCa = p["gNab"], p["gCab"], p["gNaK"], p["gNaCa"], p["fCa"]
    eta, ksat = p["eta_NaCa"], p["ksat_NaCa"]

    for i in range(n):
        V_out[i] = V
        if i == n - 1:
            break

        t = i * dt
        phase = (t - stim_offset) % period
        Istim = stim_amp if (0.0 <= phase < stim_dur and t >= stim_offset) else 0.0

        g = pb_reduced_gates(V, p)

        # pump / exchanger (treated as constant over the step)
        fNaK = 1.0 / (1.0 + 0.1245 * math.exp(-0.1 * V / RTF) +
                      0.0365 * sigma * math.exp(-V / RTF))
        INaK = gNaK * fNaK * fNaK2
        fNaCa = naca_pref / (1.0 + ksat * math.exp((eta - 1.0) * V / RTF)) * \
                (math.exp(eta * V / RTF) * p["Nai"] ** 3 * p["Cae"] -
                 math.exp((eta - 1.0) * V / RTF) * p["Nae"] ** 3 * p["Cai"])
        INaCa = gNaCa * fNaCa

        gNa_eff = gNa * m ** 3 * v ** 2
        gCa_eff = gCa * g["dinf"] * f * fCa
        gto_eff = gto * g["rinf"] * to
        gK_eff = gK * X ** 2

        aK1 = 0.1 / (1.0 + math.exp(0.06 * (V - EK - 200.0)))
        bK1 = (3.0 * math.exp(2e-4 * (V - EK + 100.0)) +
               math.exp(0.1 * (V - EK - 10.0))) / (1.0 + math.exp(-0.5 * (V - EK)))
        gK1_eff = gK1 * (aK1 / (aK1 + bK1))

        gTot = gNa_eff + gCa_eff + gto_eff + gK_eff + gK1_eff + gNab + gCab
        ERevSum = (gNa_eff * ENa + gCa_eff * ECa + gto_eff * Eto +
                   gK_eff * EK + gK1_eff * EK + gNab * ENa + gCab * ECa)

        Vinf = (ERevSum - INaK - INaCa + Istim) / gTot
        tauV = 1.0 / gTot

        m = g["minf"] + (m - g["minf"]) * math.exp(-dt / g["taum"])
        v = g["vinf"] + (v - g["vinf"]) * math.exp(-dt / g["tauv"])
        f = g["finf"] + (f - g["finf"]) * math.exp(-dt / g["tauf"])
        to = g["toinf"] + (to - g["toinf"]) * math.exp(-dt / g["tauto"])
        X = g["Xinf"] + (X - g["Xinf"]) * math.exp(-dt / g["tauX"])
        V = Vinf + (V - Vinf) * math.exp(-dt / tauV)

    t_arr = np.arange(n) * dt
    return t_arr, V_out


# -------------------------------------------------------------- AP features
def ap_features(t, V, beat_start, beat_end):
    """Features of the beat in [beat_start, beat_end) ms."""
    sel = (t >= beat_start) & (t < beat_end)
    tb, Vb = t[sel], V[sel]
    if len(Vb) < 10:
        return None

    dV = np.gradient(Vb, tb)
    i_up = int(np.argmax(dV))
    t_up = tb[i_up]
    dVdt_max = float(dV[i_up])

    V_rest = float(np.min(Vb))
    i_peak = int(np.argmax(Vb))
    V_peak = float(Vb[i_peak])
    amp = V_peak - V_rest

    # A real AP needs a fast upstroke and a decent overshoot. Without them
    # the cell failed to capture (or failed to repolarize) and every APD
    # number below is meaningless -- flag it instead of reporting it.
    captured = (dVdt_max > 50.0) and (amp > 80.0) and (V_peak > 0.0)

    out = {"V_rest": V_rest, "V_peak": V_peak, "dVdt_max": dVdt_max,
           "t_up": float(t_up), "captured": captured}
    if not captured:
        for pct in (30, 50, 90):
            out[f"APD{pct}"] = float("nan")
        return out

    for pct in (30, 50, 90):
        thresh = V_peak - (pct / 100.0) * amp
        after = np.where((tb > tb[i_peak]) & (Vb <= thresh))[0]
        out[f"APD{pct}"] = float(tb[after[0]] - t_up) if len(after) else float("nan")
    return out


if __name__ == "__main__":
    import time
    p = pb_reduced_params("epi")
    t0 = time.time()
    t, V = simulate(p, t_end=5000.0)
    print(f"simulated 5 beats in {time.time()-t0:.1f} s")
    feat = ap_features(t, V, 4050.0, 5000.0)
    for k, val in feat.items():
        print(f"  {k:10s} {val:10.2f}")
