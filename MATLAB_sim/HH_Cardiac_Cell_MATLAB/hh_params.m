function p = hh_params()
p.C     = 10;       % membrane capacitance          [nF/mm^2]
p.gL    = 0.003e3;  % leak maximal conductance      [uS/mm^2]
p.EL    = -54.387;  % leak reversal potential       [mV]
p.gK    = 0.36e3;   % potassium maximal conductance [uS/mm^2]
p.EK    = -77;      % potassium reversal potential  [mV]
p.gNa   = 1.2e3;    % sodium maximal conductance    [uS/mm^2]
p.ENa   = 50;       % sodium reversal potential     [mV]
p.Vrest = -65;      % resting (initial) potential   [mV]
end
