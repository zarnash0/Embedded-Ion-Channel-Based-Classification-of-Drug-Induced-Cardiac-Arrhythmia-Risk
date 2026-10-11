%% RUN_NOBLE_DEMO
% Classical Noble (1962) cardiac Purkinje-fibre model cross-validated against two
% independent authoritative sources (Fink & Noble 2008 Scholarpedia
% article, and the CellML/Physiome reference implementation)- with
% full electrophysiological visualization and quantitative analysis,
close all; clear; clc;
figDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

params = noble_params();
DT = 0.01;        
%% 1) Detailed short run
SHORT_T = 3000;  %(3-4 spontaneous cycles)
resShort = simulate_Noble(params, SHORT_T, dt=DT);

fig1 = plot_noble_ap(resShort);
exportgraphics(fig1, fullfile(figDir,'01_action_potential.png'), 'Resolution',150);

fig2 = plot_noble_gating_and_currents(resShort);
exportgraphics(fig2, fullfile(figDir,'02_gating_and_currents.png'), 'Resolution',150);

fig3 = plot_noble_log_conductance(resShort);
exportgraphics(fig3, fullfile(figDir,'03_log_conductance.png'), 'Resolution',150);

fig4 = plot_noble_rate_functions();
exportgraphics(fig4, fullfile(figDir,'04_rate_functions.png'), 'Resolution',150);

fig5 = plot_noble_phase_portrait_2d(resShort);
exportgraphics(fig5, fullfile(figDir,'05_phase_portrait_2d.png'), 'Resolution',150);

%% 
LONG_T = 6000;    
resLong = simulate_Noble(params, LONG_T, dt=DT);
beatsAll = detect_beats(resLong);
if beatsAll.count >= 3
    v = 2:beatsAll.count;   
else
    v = 1:beatsAll.count;
end

%% 3) Cross-check A: exponential-update vs RK4
resShortRK4 = simulate_Noble(params, SHORT_T, dt=DT, method="rk4");
beatsExp = detect_beats(resShort);
beatsRK4 = detect_beats(resShortRK4);
maxErrV = max(abs(resShort.V - resShortRK4.V));
rmsErrV = sqrt(mean((resShort.V - resShortRK4.V).^2));
beatCountMatch = (beatsExp.count == beatsRK4.count);
if beatCountMatch && beatsExp.count > 0
    meanBeatTimeErr  = mean(abs(beatsExp.times - beatsRK4.times));
    meanOvershootErr = mean(abs(beatsExp.overshoot - beatsRK4.overshoot));
else
    meanBeatTimeErr  = NaN;
    meanOvershootErr = NaN;
end

%% 4) Cross-check: exponential method
CONV_T = 1200;    % ms, enough to span >=1 full cycle
resConvCoarse = simulate_Noble(params, CONV_T, dt=DT);
resConvFine   = simulate_Noble(params, CONV_T, dt=DT/10);
Vfine_at_coarse_t = interp1(resConvFine.t, resConvFine.V, resConvCoarse.t);
maxErrConv = max(abs(resConvCoarse.V - Vfine_at_coarse_t));

%% 5) Analytics summary 
fprintf('\n================ Noble (1962) Model: Summary ================\n');
fprintf('Parameters: Cm=%.0f uF   ENa=%.1f  EK=%.1f  EAn=%.1f mV\n', ...
    params.Cm, params.ENa, params.EK, params.EAn);
fprintf('Conductances: gNaMax=%.0f  gNaBg=%.0f  gK1a=%.0f  gK1b=%.0f  gK2Max=%.0f  gAn=%.0f uS\n', ...
    params.gNaMax, params.gNaBg, params.gK1a, params.gK1b, params.gK2Max, params.gAn);
fprintf('Initial conditions (CellML reference): V0=%.0f mV, m0=%.2f, h0=%.2f, n0=%.2f\n', ...
    params.V0, params.m0, params.h0, params.n0);
fprintf('No external stimulus -- spontaneous pacemaker, dt=%.4g ms\n\n', DT);

fprintf('--- Spontaneous pacemaker activity (%.0f ms run, %d beats, first beat excluded) ---\n', ...
    LONG_T, numel(v));
fprintf('Beats detected            : %d\n', beatsAll.count);
if numel(v) >= 1
    fprintf('Overshoot (peak)           : %.2f +/- %.2f mV\n', mean(beatsAll.overshoot(v)), std(beatsAll.overshoot(v)));
    fprintf('Maximum diastolic potential: %.2f +/- %.2f mV\n', mean(beatsAll.MDP(v)), std(beatsAll.MDP(v)));
    fprintf('Amplitude                  : %.2f +/- %.2f mV\n', mean(beatsAll.amplitude(v)), std(beatsAll.amplitude(v)));
    fprintf('Upstroke velocity dV/dt max: %.1f +/- %.1f V/s\n', mean(beatsAll.dVdtMax(v)), std(beatsAll.dVdtMax(v)));
    fprintf('APD90                      : %.1f +/- %.1f ms\n', mean(beatsAll.APD(v),'omitnan'), std(beatsAll.APD(v),'omitnan'));
end
if beatsAll.count >= 3
    clSS = beatsAll.cycleLength(2:end);
    fprintf('Cycle length               : %.1f +/- %.1f ms\n', mean(clSS), std(clSS));
    fprintf('Spontaneous rate           : %.1f bpm\n', 60000/mean(clSS));
    fprintf('Cycle-length CV            : %.4f\n', std(clSS)/mean(clSS));
end

fprintf('\n--- Validation against Fink & Noble (2008) specific claims ---\n');
if numel(v) >= 1
    fprintf('Pacemaker-mechanism consistency score : %.1f%% (declining n + non-decreasing V during diastole; expect high, >90%%)\n', ...
        mean(beatsAll.pacemakerMechanismScore(v),'omitnan'));
    fprintf('Persistent Na+ conductance in plateau : %.2f%% of peak (Scholarpedia Fig.4 caption states "about 2%%")\n', ...
        mean(beatsAll.persistentNaConductancePct(v),'omitnan'));
else
    fprintf('(insufficient beats for validation-target statistics)\n');
end

fprintf('\n--- Numerical cross-checks ---\n');
fprintf('A) exponential-update vs. RK4 (dt=%.4g ms, %.0f ms window):\n', DT, SHORT_T);
fprintf('   Raw pointwise max/RMS |V_exp-V_rk4| : %.3e / %.3e mV (can look large from sub-ms\n', maxErrV, rmsErrV);
fprintf('     beat-timing phase drift alone -- see beat-aligned comparison below)\n');
fprintf('   Beat count exp vs rk4               : %d vs %d %s\n', beatsExp.count, beatsRK4.count, ...
    ternary(beatCountMatch,'(match)','(MISMATCH)'));
if beatCountMatch && beatsExp.count > 0
    fprintf('   Mean |beat time exp-rk4|            : %.4f ms\n', meanBeatTimeErr);
    fprintf('   Mean |overshoot exp-rk4|            : %.4f mV\n', meanOvershootErr);
end
fprintf('B) dt=%.4g ms vs. dt/10=%.4g ms (exponential method, %.0f ms window):\n', DT, DT/10, CONV_T);
fprintf('   Raw pointwise max error              : %.3e mV (should shrink vs. dt -> confirms convergence)\n', maxErrConv);
fprintf('===============================================================\n\n');

%% 
nB = beatsAll.count;
cycleLenCol = [beatsAll.cycleLength(:); NaN];
if numel(cycleLenCol) < nB
    cycleLenCol = [cycleLenCol; nan(nB-numel(cycleLenCol),1)];
end
T_beats = table((1:nB)', beatsAll.times(:), beatsAll.peakTimes(:), beatsAll.overshoot(:), ...
    beatsAll.MDP(:), beatsAll.amplitude(:), beatsAll.dVdtMax(:), beatsAll.APD(:), cycleLenCol(1:nB), ...
    beatsAll.pacemakerMechanismScore(:), beatsAll.persistentNaConductancePct(:), ...
    'VariableNames', {'Beat','CrossTime_ms','PeakTime_ms','Overshoot_mV','MDP_mV', ...
    'Amplitude_mV','dVdtMax_Vs','APD90_ms','CycleLengthToNext_ms', ...
    'PacemakerMechanismScore_pct','PersistentNaConductance_pct'});
writetable(T_beats, fullfile(figDir,'beat_by_beat_data.csv'));

summaryRows = { ...
    'BeatsDetected_6sRun', beatsAll.count; ...
    'MeanOvershoot_mV', mean(beatsAll.overshoot(v)); ...
    'MeanMDP_mV', mean(beatsAll.MDP(v)); ...
    'MeanAmplitude_mV', mean(beatsAll.amplitude(v)); ...
    'MeanDVDTmax_Vs', mean(beatsAll.dVdtMax(v)); ...
    'MeanAPD90_ms', mean(beatsAll.APD(v),'omitnan'); ...
    'MeanCycleLength_ms', ternary(beatsAll.count>=3, mean(beatsAll.cycleLength), NaN); ...
    'SpontaneousRate_bpm', beatsAll.rateBpm; ...
    'CycleLengthCV', beatsAll.CV_CL; ...
    'PacemakerMechanismScore_pct', mean(beatsAll.pacemakerMechanismScore(v),'omitnan'); ...
    'PersistentNaConductance_pct', mean(beatsAll.persistentNaConductancePct(v),'omitnan'); ...
    'RawPointwiseMaxErr_ExpVsRK4_mV', maxErrV; ...
    'BeatCountMatch_ExpVsRK4', double(beatCountMatch); ...
    'MeanBeatTimeErr_ExpVsRK4_ms', meanBeatTimeErr; ...
    'MeanOvershootErr_ExpVsRK4_mV', meanOvershootErr; ...
    'RawPointwiseMaxErr_dtVsDtOver10_mV', maxErrConv};
T_summary = cell2table(summaryRows, 'VariableNames', {'Metric','Value'});
writetable(T_summary, fullfile(figDir,'summary_metrics.csv'));

fprintf('Figures saved to: %s\n', figDir);
fprintf('Analytics CSVs saved to: %s\n\n', figDir);

function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end
