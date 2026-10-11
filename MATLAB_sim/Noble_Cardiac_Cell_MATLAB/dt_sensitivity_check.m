%% DT_SENSITIVITY_CHECK
% Diagnostic: checks how the measured
% peak upstroke velocity and peak effective Na+ conductance
% depend on the integration step dt
clear; clc;
params = noble_params();
T = 900; % ms, enough for one full beat plus margin

dts = [0.01, 0.005, 0.002, 0.001, 0.0005];
fprintf('%10s %14s %14s %10s\n', 'dt(ms)', 'dVdtMax(V/s)', 'peakGNa(uS)', 'nSteps');
for i = 1:numel(dts)
    dt = dts(i);
    res = simulate_Noble(params, T, dt=dt);
    beats = detect_beats(res);
    if beats.count >= 2
        dvdt = beats.dVdtMax(2);
        iFrom = find(res.t >= beats.times(2), 1, 'first');
        peakGNa = max(res.gNa(iFrom:end));
    else
        dvdt = NaN;
        peakGNa = NaN;
    end
    fprintf('%10.4f %14.2f %14.1f %10d\n', dt, dvdt, peakGNa, numel(res.t));
end
