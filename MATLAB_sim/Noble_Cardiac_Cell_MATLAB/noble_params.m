function p = noble_params()

%   Equations/constants cross-validated against two independent sources:
%     (1) Fink M & Noble D (2008). "Noble model." Scholarpedia 3(2):1803.
%         doi:10.4249/scholarpedia.1803 -- states "Parameters are given
%         in millivolt, microsiemens and seconds."
%     (2) The CellML/Physiome Model Repository reference implementation
%         (models.physiomeproject.org/workspace/noble_1962), built
%         directly from the 1962 paper 

p.Cm  = 12;      % membrane capacitance                          [uF]
p.ENa = 40;      % sodium reversal potential                     [mV]
p.EK  = -100;    % potassium reversal potential (shared iK1/iK2)  [mV]
p.EAn = -60;     % anion (chloride) leak reversal potential       [mV]

p.gNaMax  = 400000;  % fast, HH-style gated Na+ conductance coeff. (m^3 h) [uS]
p.gNaBg   = 140;     % persistent/background Na+ conductance (non-gated)   [uS]
p.gK1a    = 1200;    % i_K1 (anomalous/inward rectifier) coefficient 1     [uS]
p.gK1b    = 15;      % i_K1 (anomalous/inward rectifier) coefficient 2     [uS]
p.gK2Max  = 1200;    % i_K2 (slow delayed-rectifier-like, gated by n^4)    [uS]
p.gAn     = 75;      % anion/chloride leak conductance (ohmic)             [uS]

% Reference initial conditions (CellML/Physiome), NOT recomputed steady
% states -- used verbatim to match the authoritative reference run.
p.V0 = -87;    % mV
p.m0 = 0.01;   
p.h0 = 0.8;    
p.n0 = 0.01;   
end
