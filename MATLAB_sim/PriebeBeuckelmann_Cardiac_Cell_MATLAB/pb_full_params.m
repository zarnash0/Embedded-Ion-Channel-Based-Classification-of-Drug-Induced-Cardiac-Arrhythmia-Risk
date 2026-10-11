function p = pb_full_params(condition)

if nargin < 1, condition = "nonfailing"; end
condition = string(condition);
p.condition = condition;
failing = (condition == "failing");

p.Cm = 153.4;    % pF, total cell membrane capacitance
cst = pb_constants();
p.R = cst.R;  p.T = cst.T;  p.F = cst.F;  p.RTF = cst.RTF;

p.gNa    = 16.0;      % mS/uF
p.gCa    = 0.064;
p.gKs    = 0.02;
p.gKr    = 0.015;
if failing
    p.gto  = 0.191;  p.gK1 = 2.0;  p.gNab = 0;      p.gCab = 0.0013;
    p.INaKmax = 0.75; p.KNaCa = 1650;
else
    p.gto  = 0.3;    p.gK1 = 2.5;  p.gNab = 0.001;  p.gCab = 0.00085;
    p.INaKmax = 1.3;  p.KNaCa = 1000;
end

p.KmNai_NaK = 10;   p.KmKo_NaK = 1.5;
p.KmNa_NaCa = 87.5; p.KmCa_NaCa = 1.38; p.ksat_NaCa = 0.1; p.eta_NaCa = 0.35;
p.KmCa_fCa  = 0.0006;   % ICa Ca2+-dependent inactivation factor fCa

% Extracellular / fixed concentrations (mM)
p.Nao = 138;  p.Ko = 4;  p.Cao = 2;

% CICR (Ca2+-induced Ca2+ release)- the "normal, AP-triggered" pathway
% only (the CellML reference also contains a "spontaneous release /
% overload" pathway used solely for Priebe & Beuckelmann's separate DAD
% experiments in their Figs. 14-15, wired through 3 additional tracking
% ODEs (OVRLDtrack, OVRLDtrack2, OVRLDtrack3) that are multiplied by
% 0 everywhere they enter the active computation in the baseline
% reference encoding- i.e. that pathway is inert for ordinary paced APs
% and is omitted 
p.Grel_max = 22;     p.Kmrel = 0.0008;   p.delta_Ca_ith = 5e-6;

if failing
    p.Iupmax = 0.0015; p.Kleak = 0.00017;
else
    p.Iupmax = 0.0045; p.Kleak = 0.00026;
end
p.Kmup  = 0.00092;
p.tautr = 180;

p.KmTn = 0.0005;  p.KmCMDN = 0.00238;  p.Tnmax = 0.07;  p.CMDNmax = 0.05;  p.buffon = 1;
p.CSQNmax = 10;   p.KmCSQN = 0.8;

p.preplength = 100;  p.radius = 11;
p.volume = pi*p.preplength*p.radius^2;
p.Vmyo = 0.68*p.volume;  p.VJSR = 0.0048*p.volume;  p.VNSR = 0.0552*p.volume;

p.EtoNum = 0.043*p.Nao + p.Ko;          % Eto numerator
p.EKsNum = 0.01833*p.Nao + p.Ko;        % EKs numerator
p.dCoefA = 14.9859/(16.6813*sqrt(2*pi));
p.dCoefB = 5.3/(14.93*sqrt(2*pi));
p.sigmaNaK0365 = 0.0365*(1/7)*(exp(p.Nao/67.3)-1);   % 0.0365*sigma
p.KoFrac  = p.Ko/(p.Ko+p.KmKo_NaK);
p.Nao3    = p.Nao^3;
p.NaCaDen = (p.KmNa_NaCa^3 + p.Nao^3)*(p.KmCa_NaCa + p.Cao);

% Reference initial conditions (CellML, nonfailing baseline; used
p.V0 = -90.7796417483135;
p.m0 = 0.000585525582501575;
p.h0 = 0.995865529216237;
p.j0 = 0.997011204496203;
p.d0 = 2.50653215966786e-10;
p.f0 = 0.92130376850548;
p.r0 = 1.75032478501027e-5;
p.to0 = 0.999897251531651;
p.Xs0 = 0.00885658064818147;
p.Xr0 = 0.000215523048438941;
p.Nai0 = 10;   p.Ki0 = 140;
% SR calcium load is CONDITION-SPECIFIC. Priebe & Beuckelmann state it
% directly: simulations "were started with [Ca2+]NSR = [Ca2+]JSR = 2.5
% mmol/L in a nonfailing myocyte and [Ca2+]NSR = [Ca2+]JSR = 1.0 mmol/L
% in a failing myocyte on the basis of experimental data" - the failing
% myocyte's SR is depleted (reduced SERCA activity; the paper's Iupmax
% drops 0.0045 -> 0.0015 for the same reason).
p.Cai0 = 0.0002;
if failing
    p.CaJSR0 = 1.0;  p.CaNSR0 = 1.0;
else
    p.CaJSR0 = 2.5;  p.CaNSR0 = 2.5;
end
p.APtrack0 = 0;  p.APtrack20 = 0;  p.Cainfluxtrack0 = 0;
end
