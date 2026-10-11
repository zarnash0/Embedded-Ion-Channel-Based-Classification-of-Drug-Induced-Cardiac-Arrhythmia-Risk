function pal = noble_palette()
pal.V    = [0.10 0.10 0.10];   % membrane potential   - near black
pal.stim = [0.55 0.55 0.55];   % stimulus / perturbation current - grey
pal.m    = [0.00 0.53 0.30];   % Na+ activation   (m)  - green
pal.h    = [0.80 0.15 0.15];   % Na+ inactivation (h)  - red
pal.n    = [0.00 0.45 0.74];   % K+  activation   (n)  - blue
pal.iNa  = [0.00 0.53 0.30];   % sodium current / conductance
pal.iK1  = [0.49 0.18 0.56];   % anomalous (inward) rectifier - purple
pal.iK2  = [0.00 0.45 0.74];   % slow delayed-rectifier-like K+ current
pal.iAn  = [0.60 0.60 0.60];   % anion/chloride leak
pal.depol   = [0.00 0.45 0.85];  % dV/dt > 0 (depolarizing)  - blue, Fenton/Virtual-Heart convention
pal.repol   = [0.86 0.20 0.55];  % dV/dt < 0 (repolarizing)  - pink, Fenton/Virtual-Heart convention
pal.grid = [0.85 0.85 0.85];
end
