function p = pb_reduced_params(cellType)
if nargin < 1, cellType = "epi"; end
cellType = string(cellType);

p.cellType = cellType;

p.gNa    = 16.0;     % nS/pF (paper's units; current is g*gates*(V-E) in pA/pF)
p.gCa    = 0.064;
p.gNab   = 0.001;
p.gCab   = 0.00085;
p.gNaK   = 1.3;       % pA/pF scaling factor (Bernus Table 1 labels it nS/pF, but
                      % it is used exactly as PB's I-bar_NaK, a current-density scale)
p.gNaCa  = 1000;      % likewise, PB's kNaCa scaling factor (pA/pF-equivalent)

switch cellType
    case "epi"
        p.gto = 0.4; p.gK = 0.019; p.invP = 1;   p.Vshift = 0;
    case "M"
        p.gto = 0.35; p.gK = 0.013; p.invP = 1.7; p.Vshift = -4;
    case "endo"
        p.gto = 0.13; p.gK = 0.019; p.invP = 2.8; p.Vshift = -12;
    otherwise
        error('pb_reduced_params:cellType', 'cellType must be "epi", "M", or "endo".');
end
p.gK1 = 3.9;

p.Nai = 10;      p.Nae = 138;    % mM
p.Ki  = 140;     p.Ke  = 4;      % mM
p.Cai = 0.0004;  p.Cae = 2;      % mM  (0.0004 mM = 400 nM)
p.fCa = 0.6;     % = 1/(1+Cai/0.0006) evaluated at the fixed Cai above

p.KmNai_NaK = 10;    p.KmKo_NaK = 1.5;             % mM, INaK
p.KmNa_NaCa = 87.5;  p.KmCa_NaCa = 1.38;           % mM, INaCa
p.ksat_NaCa = 0.1;   p.eta_NaCa  = 0.35;

p.V0 = -90.2;    % mV, resting potential reported for the reduced model (Table 3)
% Gating variables are not hardcoded here- SIMULATE_PB_REDUCED
% initializes them at their steady-state value for V0, computed from the
% actual gate equations (PB_REDUCED_GATES), as done for the
% earlier Hodgkin-Huxley and Noble models in this project.
end
