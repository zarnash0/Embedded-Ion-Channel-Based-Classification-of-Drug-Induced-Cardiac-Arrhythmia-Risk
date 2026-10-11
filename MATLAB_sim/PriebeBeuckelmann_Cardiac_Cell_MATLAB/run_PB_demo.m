%% RUN_PB_DEMO
% Priebe & Beuckelmann (1998) human ventricular myocyte model and the
% Bernus et al. (2002) six-variable reduction of it -- the "third part" of
% the hardware cardiac-cell project (Hodgkin-Huxley -> Noble 1962 ->
% Priebe-Beuckelmann). This script reproduces, in order: the full model's
% nonfailing action potential with its currents, gates and Ca2+ handling;
% the paper's central heart-failure result (AP prolongation in the failing
% myocyte); Bernus's transmural cell types (epi / M / endo); a full-vs-
% reduced model comparison; and APD restitution for both models.
%
%   Priebe L & Beuckelmann DJ (1998). Circ Res 82:1206-1223.
%   Bernus O, Wilders R, Zemlin CW, Verschelde H, Panfilov AV (2002).
%       Am J Physiol Heart Circ Physiol 282:H2296-H2308.
%
% RUNTIME. The full model is RK4 at dt = 0.005 ms and costs roughly 30 s
% of wall clock per 1200 ms of simulated time, i.e. ~2 min per 5-s run.
% Sections 1, 2, 4 and 5 each pay that price, so budget on the order of
% 15 minutes for the whole script; the reduced model is essentially free
% by comparison (that contrast is itself one of the results). Progress
% messages are printed before every long block, and after the first run
% the script prints a measured projection of the total remaining time.
%
% Companion plotting functions used below (single-run figures):
%   PLOT_PB_AP, PLOT_PB_CURRENTS, PLOT_PB_GATES, PLOT_PB_CALCIUM,
%   PLOT_PB_PHASE3D, PLOT_PB_COMPARISON.
% The three multi-run overlays (failing vs nonfailing, the three cell
% types, and the restitution curves) are local functions at the bottom of
% this file, since they compose several simulations and exist only for
% this demo.

close all; clear; clc;
figDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

pal = pb_palette();

DT_FULL = 0.005;   % ms, full model  (its own validated default step)
DT_RED  = 0.02;    % ms, reduced model (its own validated default step)
T_MAIN  = 5000;    % ms, 5 paced beats at 1 Hz; the LAST beat is analysed
T_CMP   = 3000;    % ms, shorter window for the full-vs-reduced timing run

% Restitution protocol for the full model (Section 5). Declared here so
% Section 1 can project the total runtime from its own measured cost.
NS1_FULL = 3;                          % S1 beats (10 would take over an hour)
DI_FULL  = [30 60 120 250 400];        % ms, diastolic intervals
DI_RED   = [20 30 50 75 100 150 200 300 400 500];

fprintf('\n=========== Priebe-Beuckelmann (1998) / Bernus (2002) demo ===========\n');
fprintf('Full model : RK4, dt = %.4g ms, %.0f ms of model time per run\n', DT_FULL, T_MAIN);
fprintf('Reduced    : Rush-Larsen, dt = %.4g ms, %.0f ms per run\n', DT_RED, T_MAIN);
fprintf('Pacing     : 1 Hz, 3 ms pulses of 15 pA/pF (both papers'' protocol)\n');
fprintf('This takes several minutes -- the full model is the expensive part.\n\n');

%% 0) JIT warm-up ---------------------------------------------------------
% Without this the first timed call below is the first simulate_PB_full of
% the session and is charged for MATLAB's just-in-time compilation as well
% as the arithmetic. Those runtimes feed speedupWall and are archived to
% pb_efficiency_metrics.csv, so an unwarmed first pass silently inflates
% the full model's cost and the reported speedup. COMPARE_ALL_MODELS was
% already fixed for exactly this; RUN_PB_DEMO had not been.
fprintf('[0/6] Warming up the JIT (untimed)...\n');
wf = simulate_PB_full(pb_full_params("nonfailing"), 5, dt=DT_FULL);      %#ok<NASGU>
wr = simulate_PB_reduced(pb_reduced_params("epi"), 5, dt=DT_RED);        %#ok<NASGU>
clear wf wr

%% 1) Full model, nonfailing ---------------------------------------------
fprintf('[1/6] Full model, NONFAILING: simulating %.0f ms at dt = %.4g ms...\n', T_MAIN, DT_FULL);
paramsNF = pb_full_params("nonfailing");
tic;
resNF = simulate_PB_full(paramsNF, T_MAIN, dt=DT_FULL);
runtimeFullNF = toc;
apNF = detect_pb_ap(resNF);
mNF  = localLastBeatMetrics(resNF, apNF);
fprintf('      done in %.1f s (%.4f s per ms of model time).\n', ...
    runtimeFullNF, runtimeFullNF/T_MAIN);

% Project the rest of the run from the cost just measured. The restitution
% term counts the S1-only reference run plus one run per DI; 400 ms stands
% in for the last S1 beat's APD90, which sets where the S2 lands. Assumes
% COMPARE_PB_MODELS integrates the full model once over T_CMP.
restFullModelMs = (50 + (NS1_FULL-1)*1000 + 1000) ...
                + sum(50 + (NS1_FULL-1)*1000 + 400 + DI_FULL + 800);
projModelMs = T_MAIN + T_CMP + restFullModelMs;
fprintf('      remaining full-model work: %.1f s of model time -> ~%.0f min projected.\n', ...
    projModelMs/1000, runtimeFullNF/T_MAIN*projModelMs/60);

fig1 = plot_pb_ap(resNF);
exportgraphics(fig1, fullfile(figDir,'01_full_nonfailing_action_potential.png'), 'Resolution',150);

fig2 = plot_pb_currents(resNF);
exportgraphics(fig2, fullfile(figDir,'02_full_nonfailing_currents.png'), 'Resolution',150);

fig3 = plot_pb_gates(resNF);
exportgraphics(fig3, fullfile(figDir,'03_full_nonfailing_gates.png'), 'Resolution',150);

fig4 = plot_pb_calcium(resNF);
exportgraphics(fig4, fullfile(figDir,'04_full_nonfailing_calcium.png'), 'Resolution',150);

fig5 = plot_pb_phase3d(resNF);
exportgraphics(fig5, fullfile(figDir,'05_full_nonfailing_phase_portrait_3d.png'), 'Resolution',150);

%% 2) Full model, FAILING vs nonfailing ----------------------------------
% This is the paper's central result: in the failing myocyte the AP is
% markedly prolonged, driven by enhanced INaCa, slowed [Ca2+]i decay, and
% reduced IK1 and INaK (all encoded in PB_FULL_PARAMS("failing")).
fprintf('[2/6] Full model, FAILING: simulating %.0f ms at dt = %.4g ms (~%.0f s)...\n', ...
    T_MAIN, DT_FULL, runtimeFullNF);
paramsF = pb_full_params("failing");
tic;
resF = simulate_PB_full(paramsF, T_MAIN, dt=DT_FULL);
runtimeFullF = toc;
apF = detect_pb_ap(resF);
mF  = localLastBeatMetrics(resF, apF);
fprintf('      done in %.1f s.\n', runtimeFullF);

fig6 = localPlotFailingOverlay(resNF, resF, apNF, apF, pal);
exportgraphics(fig6, fullfile(figDir,'06_failing_vs_nonfailing.png'), 'Resolution',150);

% Published reference values (Priebe & Beuckelmann 1998, 1 Hz steady state)
refNF = struct('RMP',-89.7, 'plateau',19.0, 'dVdtMax',NaN, ...
               'APD25',262.9, 'APD50',310.2, 'APD90',374.0, ...
               'CaiRest',120, 'CaiPeak',614);
refF  = struct('RMP',-85.6, 'plateau',29.2, 'dVdtMax',NaN, ...
               'APD25',305.7, 'APD50',374.5, 'APD90',548.8, ...
               'CaiRest',NaN, 'CaiPeak',334);

fprintf('\n      --- Heart-failure result: this implementation vs. the paper ---\n');
fprintf('      %-12s %10s %10s %10s %10s\n', '', 'APD25', 'APD50', 'APD90', 'RMP');
fprintf('      %-12s %10.1f %10.1f %10.1f %10.1f\n', 'nonfailing', mNF.APD25, mNF.APD50, mNF.APD90, mNF.RMP);
fprintf('      %-12s %10.1f %10.1f %10.1f %10.1f  <- Priebe & Beuckelmann 1998\n', ...
    ' (published)', refNF.APD25, refNF.APD50, refNF.APD90, refNF.RMP);
fprintf('      %-12s %10.1f %10.1f %10.1f %10.1f\n', 'failing', mF.APD25, mF.APD50, mF.APD90, mF.RMP);
fprintf('      %-12s %10.1f %10.1f %10.1f %10.1f  <- Priebe & Beuckelmann 1998\n', ...
    ' (published)', refF.APD25, refF.APD50, refF.APD90, refF.RMP);
fprintf('      APD90 prolongation in failure: %+.1f%% here, %+.1f%% published.\n\n', ...
    100*(mF.APD90/mNF.APD90 - 1), 100*(refF.APD90/refNF.APD90 - 1));

%% 3) Reduced model, three cell types ------------------------------------
% Mirrors Bernus et al. Fig. 6A-C: epicardial, midmyocardial (M) and
% endocardial cells differ only in gto, gK and the Ito kinetics.
fprintf('[3/6] Reduced model: simulating epi, M and endo (%.0f ms each, dt = %.4g ms)...\n', ...
    T_MAIN, DT_RED);
paramsEpi  = pb_reduced_params("epi");
paramsM    = pb_reduced_params("M");
paramsEndo = pb_reduced_params("endo");

tic;
resEpi = simulate_PB_reduced(paramsEpi, T_MAIN, dt=DT_RED);
runtimeRedEpi = toc;                      % timed for the efficiency table
resM    = simulate_PB_reduced(paramsM,    T_MAIN, dt=DT_RED);
resEndo = simulate_PB_reduced(paramsEndo, T_MAIN, dt=DT_RED);

apEpi  = detect_pb_ap(resEpi);   mEpi  = localLastBeatMetrics(resEpi,  apEpi);
apM    = detect_pb_ap(resM);     mM    = localLastBeatMetrics(resM,    apM);
apEndo = detect_pb_ap(resEndo);  mEndo = localLastBeatMetrics(resEndo, apEndo);
fprintf('      done (epi run: %.2f s, i.e. %.0fx faster than the full model''s %.1f s).\n', ...
    runtimeRedEpi, runtimeFullNF/runtimeRedEpi, runtimeFullNF);

fig7 = localPlotCellTypes({resEpi, resM, resEndo}, {mEpi, mM, mEndo}, ...
                          {apEpi, apM, apEndo}, {'Epicardial','Midmyocardial (M)','Endocardial'}, pal);
exportgraphics(fig7, fullfile(figDir,'07_reduced_cell_types.png'), 'Resolution',150);

%% 4) Full vs reduced comparison -----------------------------------------
% Shorter window here (T_CMP) purely to keep the wall clock sane: this
% block pays for a second full-model run on top of Sections 1 and 2.
fprintf('[4/6] Full vs reduced comparison over %.0f ms (~%.0f s, dominated by the full model)...\n', ...
    T_CMP, runtimeFullNF*T_CMP/T_MAIN);
cmp = compare_pb_models(tEnd=T_CMP, dtFull=DT_FULL, dtReduced=DT_RED, ...
                        condition="nonfailing", cellType="epi");
fig8 = plot_pb_comparison(cmp);
exportgraphics(fig8, fullfile(figDir,'08_full_vs_reduced.png'), 'Resolution',150);

%% 5) APD restitution, both models ---------------------------------------
% The reduced model is cheap enough for the textbook protocol (10 S1
% beats, 10 diastolic intervals). The full model is not: 10 S1 beats x 10
% DIs would be well over an hour, so it gets a deliberately modest
% protocol. That is an honest limitation, not a modelling choice -- and
% note Bernus's own observation that the full model needs >1 minute of
% model time to re-equilibrate after a rate change, so a 3-beat S1 train
% leaves it further from steady state than the reduced model's 10-beat
% train leaves that one. The two curves are therefore comparable in shape
% but not perfectly matched in conditioning.
fprintf('[5/6] APD restitution, reduced model (10 S1 beats, %d diastolic intervals)...\n', numel(DI_RED));
restRed = compute_pb_restitution(@simulate_PB_reduced, paramsEpi, ...
    nS1=10, S1period=1000, DIlist=DI_RED, dt=DT_RED);

fprintf('[5/6] APD restitution, full model (%d S1 beats, %d diastolic intervals) -- the slow one...\n', ...
    NS1_FULL, numel(DI_FULL));
restFull = compute_pb_restitution(@simulate_PB_full, paramsNF, ...
    nS1=NS1_FULL, S1period=1000, DIlist=DI_FULL, dt=DT_FULL);

fig9 = localPlotRestitution(restFull, restRed, pal);
exportgraphics(fig9, fullfile(figDir,'09_apd_restitution.png'), 'Resolution',150);

%% 6) Console summary and CSV export -------------------------------------
fprintf('[6/6] Writing summary...\n');

% Published reference values for the reduced model (Bernus et al. Table 3
% quotes a single set of steady-state values; we line it up against the
% default epicardial parameters, and leave M/endo blank rather than
% attribute numbers to them that the table does not separate out).
refRed = struct('RMP',-90.2, 'plateau',NaN, 'dVdtMax',386.2, ...
                'APD25',NaN, 'APD50',254.0, 'APD90',357.2, ...
                'CaiRest',400, 'CaiPeak',400);
refNone = struct('RMP',NaN, 'plateau',NaN, 'dVdtMax',NaN, ...
                 'APD25',NaN, 'APD50',NaN, 'APD90',NaN, ...
                 'CaiRest',NaN, 'CaiPeak',NaN);

hdr = {'RMP','Peak','Plateau','dV/dt','APD25','APD50','APD90','Cai_rest','Cai_peak'};
uni = {'(mV)','(mV)','(mV)','(V/s)','(ms)','(ms)','(ms)','(nM)','(nM)'};

fprintf('\n================== Priebe-Beuckelmann / Bernus: Summary ==================\n');
fprintf('Metrics are from the LAST paced beat of a %.0f ms, 1 Hz run (closest to\n', T_MAIN);
fprintf('steady state). "Plateau" is our own operational measure -- the peak of the\n');
fprintf('dome, max V over 50-150 ms after the stimulus -- not a published definition.\n\n');
localPrintHeader(hdr, uni);
localPrintRow('Full, nonfailing',   localRowVector(mNF));
localPrintRow('  published (PB98)', localRowVector(refNF));
localPrintRow('Full, failing',      localRowVector(mF));
localPrintRow('  published (PB98)', localRowVector(refF));
localPrintRow('Reduced, epi',       localRowVector(mEpi));
localPrintRow('  published (B02)',  localRowVector(refRed));
localPrintRow('Reduced, M',         localRowVector(mM));
localPrintRow('Reduced, endo',      localRowVector(mEndo));

fprintf('\n--- Heart failure (full model) ---\n');
fprintf('APD90 prolongation      : %+.1f%% (published %+.1f%%)\n', ...
    100*(mF.APD90/mNF.APD90-1), 100*(refF.APD90/refNF.APD90-1));
fprintf('Peak [Ca^2+]_i change   : %+.1f%% (published %+.1f%%)\n', ...
    100*(mF.CaiPeak/mNF.CaiPeak-1), 100*(refF.CaiPeak/refNF.CaiPeak-1));
fprintf('RMP depolarizes by      : %+.1f mV (published %+.1f mV)\n', ...
    mF.RMP-mNF.RMP, refF.RMP-refNF.RMP);
fprintf('Mechanism (per the paper): enhanced INaCa, slowed [Ca2+]i decay,\n');
fprintf('reduced IK1 and INaK -- all present in PB_FULL_PARAMS("failing").\n');

% --- Efficiency. Measured wall-clock ratio conflates two separate gains:
% fewer/cheaper equations AND a 4x larger stable time step. Reporting the
% per-step cost ratio alongside it separates them honestly; Bernus's own
% 4.9x figure is a like-for-like comparison and should not be expected to
% equal our raw wall-clock ratio.
nStepsFull = numel(resNF.t);
nStepsRed  = numel(resEpi.t);
speedupWall    = runtimeFullNF/runtimeRedEpi;
speedupPerStep = (runtimeFullNF/nStepsFull)/(runtimeRedEpi/nStepsRed);

fprintf('\n--- Computational cost (this implementation, %.0f ms of model time) ---\n', T_MAIN);
fprintf('Full model   : %8.2f s   (%d steps at dt = %.4g ms, 18 ODEs)\n', ...
    runtimeFullNF, nStepsFull, DT_FULL);
fprintf('Reduced model: %8.2f s   (%d steps at dt = %.4g ms,  6 ODEs)\n', ...
    runtimeRedEpi, nStepsRed, DT_RED);
fprintf('Wall-clock speedup      : %.1fx  (includes the 4x larger stable step)\n', speedupWall);
fprintf('Per-step cost ratio     : %.1fx  (equations only, step size factored out)\n', speedupPerStep);
fprintf('State variables         : 15 model states (18 integrated here, incl. 3\n');
fprintf('                          CICR-tracking variables) -> 6, i.e. Bernus''s 2.5x\n');
fprintf('                          memory reduction.\n');
fprintf('Bernus et al. report    : 4.9x faster than the full model, 3.8x faster than\n');
fprintf('                          the full model with [Na+]i/[K+]i clamped, memory\n');
fprintf('                          down 2.5x, no long-term ionic drift, and realistic\n');
fprintf('                          rate adaptation within a few beats (the full model\n');
fprintf('                          needs >1 min of model time to re-equilibrate after\n');
fprintf('                          a rate change; over 1000 s at 1 Hz its mean [Na+]i\n');
fprintf('                          falls 10 -> 8.1 mM and APD creeps from ~360 to\n');
fprintf('                          425 ms, never reaching steady state).\n');
fprintf('Drift check here        : [Na+]i %.2f -> %.2f mM over %.0f ms (the drift is\n', ...
    resNF.Nai(1), resNF.Nai(end), T_MAIN);
fprintf('                          slow, so a 5 s window only shows its beginning).\n');

fprintf('\n--- APD restitution (S1-S2, APD90 of the S2 beat) ---\n');
fprintf('%-28s %10s %10s\n', '', 'shortest DI', 'longest DI');
localPrintRestSummary(restFull);
localPrintRestSummary(restRed);
fprintf('Restitution is the standard validation for these models; Bernus reports his\n');
fprintf('reduced model matching the Morgan et al. (1992) human right-ventricle data\n');
fprintf('slightly better than the full model does.\n');
fprintf('=========================================================================\n\n');

% --- CSV 1: AP metrics for every configuration run ----------------------
metricRows = [ ...
    localMetricRow('PB_full',    'nonfailing', mNF,   refNF,   'Priebe & Beuckelmann 1998'); ...
    localMetricRow('PB_full',    'failing',    mF,    refF,    'Priebe & Beuckelmann 1998'); ...
    localMetricRow('PB_reduced', 'epi',        mEpi,  refRed,  'Bernus et al. 2002, Table 3'); ...
    localMetricRow('PB_reduced', 'M',          mM,    refNone, ''); ...
    localMetricRow('PB_reduced', 'endo',       mEndo, refNone, '')];
T_metrics = cell2table(metricRows, 'VariableNames', ...
    {'Model','Configuration','RMP_mV','Overshoot_mV','Amplitude_mV','Plateau_mV', ...
     'dVdtMax_Vs','APD25_ms','APD50_ms','APD90_ms','CaiRest_nM','CaiPeak_nM', ...
     'Ref_RMP_mV','Ref_Plateau_mV','Ref_dVdtMax_Vs','Ref_APD25_ms','Ref_APD50_ms', ...
     'Ref_APD90_ms','Ref_CaiRest_nM','Ref_CaiPeak_nM','Ref_Source'});
writetable(T_metrics, fullfile(figDir,'ap_metrics_all_configurations.csv'));

% --- CSV 2: comparison / efficiency numbers -----------------------------
effRows = { ...
    'ModelTime_ms',                       T_MAIN; ...
    'dt_full_ms',                         DT_FULL; ...
    'dt_reduced_ms',                      DT_RED; ...
    'Steps_full',                         nStepsFull; ...
    'Steps_reduced',                      nStepsRed; ...
    'Runtime_full_s',                     runtimeFullNF; ...
    'Runtime_full_failing_s',             runtimeFullF; ...
    'Runtime_reduced_s',                  runtimeRedEpi; ...
    'WallClockSpeedup_x',                 speedupWall; ...
    'PerStepCostRatio_x',                 speedupPerStep; ...
    'StateVariables_full_model',          15; ...
    'StateVariables_full_integrated',     18; ...
    'StateVariables_reduced',             6; ...
    'MemoryReduction_x_published',        2.5; ...
    'Bernus_Speedup_x_published',         4.9; ...
    'Bernus_Speedup_x_NaKclamped',        3.8; ...
    'dAPD90_reduced_minus_full_ms',       mEpi.APD90 - mNF.APD90; ...
    'dRMP_reduced_minus_full_mV',         mEpi.RMP   - mNF.RMP; ...
    'ddVdtMax_reduced_minus_full_Vs',     mEpi.dVdtMax - mNF.dVdtMax; ...
    'Nai_full_start_mM',                  resNF.Nai(1); ...
    'Nai_full_end_mM',                    resNF.Nai(end); ...
    'Restitution_full_nS1',               restFull.nS1; ...
    'Restitution_reduced_nS1',            restRed.nS1; ...
    'Restitution_full_APD90_shortestDI_ms',  localFirstFinite(restFull.APD90); ...
    'Restitution_full_APD90_longestDI_ms',   localLastFinite(restFull.APD90); ...
    'Restitution_reduced_APD90_shortestDI_ms', localFirstFinite(restRed.APD90); ...
    'Restitution_reduced_APD90_longestDI_ms',  localLastFinite(restRed.APD90)};

% COMPARE_PB_MODELS returns its own quantitative summary; append every
% scalar field of it generically so nothing it computes is lost here, and
% so this script does not have to track that struct's exact contents.
cmpFields = fieldnames(cmp);
for iField = 1:numel(cmpFields)
    valField = cmp.(cmpFields{iField});
    if (isnumeric(valField) || islogical(valField)) && isscalar(valField)
        effRows(end+1,:) = {['compare_' cmpFields{iField}], double(valField)}; %#ok<SAGROW>
    end
end
T_eff = cell2table(effRows, 'VariableNames', {'Metric','Value'});
writetable(T_eff, fullfile(figDir,'comparison_efficiency_metrics.csv'));

% --- CSV 3: the restitution curves themselves ---------------------------
% Not one of the two summary tables above, but the raw curve is the whole
% point of Section 5 and would otherwise survive only as a PNG.
restModel = [repmat({restFull.modelName}, numel(restFull.DI), 1); ...
             repmat({restRed.modelName},  numel(restRed.DI),  1)];
T_rest = table(restModel, ...
    [restFull.DI(:); restRed.DI(:)], ...
    [restFull.APD90(:); restRed.APD90(:)], ...
    [restFull.S2overshoot(:); restRed.S2overshoot(:)], ...
    double([restFull.captured(:); restRed.captured(:)]), ...
    [repmat(restFull.S1_APD90, numel(restFull.DI), 1); repmat(restRed.S1_APD90, numel(restRed.DI), 1)], ...
    [repmat(restFull.nS1, numel(restFull.DI), 1); repmat(restRed.nS1, numel(restRed.DI), 1)], ...
    'VariableNames', {'Model','DI_ms','S2_APD90_ms','S2_Overshoot_mV','Captured','S1_APD90_ms','nS1'});
writetable(T_rest, fullfile(figDir,'apd_restitution.csv'));

fprintf('Figures saved to: %s\n', figDir);
fprintf('  01-05 full model (nonfailing), 06 failing vs nonfailing,\n');
fprintf('  07 reduced cell types, 08 full vs reduced, 09 APD restitution\n');
fprintf('Analytics CSVs saved to: %s\n', figDir);
fprintf('  ap_metrics_all_configurations.csv, comparison_efficiency_metrics.csv,\n');
fprintf('  apd_restitution.csv\n\n');

%% ---------------------------------------------------------------- local
function m = localLastBeatMetrics(res, ap)
%LOCALLASTBEATMETRICS Metrics of the last paced beat (closest to steady
%   state) of a simulation, combining DETECT_PB_AP's output with the two
%   quantities it does not measure: the plateau potential and [Ca2+]i.
% Last beat that is actually MEASURABLE, not merely the last detected
% one: the final beat's window can be clipped by the end of the record,
% leaving NaN APDs that would otherwise propagate silently into the
% console table, both CSVs and the figure annotations -- and would
% disagree with plot_pb_ap, which independently picks the last good beat.
k = find(all(isfinite(ap.APD),2)' & ap.captured, 1, 'last');
if isempty(k)
    warning('run_PB_demo:noCompleteBeat', ...
        'No fully-measurable beat; falling back to the last detected beat (metrics may be NaN).');
    k = ap.count;
end
m.beatIndex = k;
m.RMP       = ap.RMP(k);
m.overshoot = ap.overshoot(k);
m.amplitude = ap.amplitude(k);
m.dVdtMax   = ap.dVdtMax(k);
m.APD25     = ap.APD(k,1);
m.APD50     = ap.APD(k,2);
m.APD90     = ap.APD(k,3);

% Plateau: our own operational definition -- the top of the dome, taken
% as max V over 50-150 ms after the stimulus. Neither paper states how it
% measured its quoted plateau potential, so this is flagged as ours; the
% window starts after 50 ms to skip the Ito-driven notch.
tStim = ap.stimTimes(k);
win = res.t >= tStim+50 & res.t <= tStim+150;
if any(win)
    m.plateau = max(res.V(win));
else
    m.plateau = NaN;   % record ends before the plateau window
end

if isfield(res, 'Cai')
    i0 = find(res.t >= tStim, 1, 'first');
    m.CaiRest = res.Cai(i0)*1e6;              % mM -> nM
    m.CaiPeak = max(res.Cai(i0:end))*1e6;
else
    % The reduced model has no Ca2+ state at all: a fixed 400 nM [Ca2+]i
    % is its defining simplification, so "rest" and "peak" coincide.
    m.CaiRest = res.params.Cai*1e6;
    m.CaiPeak = m.CaiRest;
end
end

function v = localRowVector(m)
%LOCALROWVECTOR Pack a metrics (or published-reference) struct into the
%   fixed column order used by the console summary table.
v = [m.RMP, localField(m,'overshoot'), m.plateau, m.dVdtMax, ...
     m.APD25, m.APD50, m.APD90, m.CaiRest, m.CaiPeak];
end

function x = localField(s, name)
if isfield(s, name), x = s.(name); else, x = NaN; end
end

function localPrintHeader(hdr, uni)
fprintf('%-20s', 'Configuration');
fprintf('%10s', hdr{:}); fprintf('\n');
fprintf('%-20s', '');
fprintf('%10s', uni{:}); fprintf('\n');
end

function localPrintRow(label, v)
%LOCALPRINTROW Print one summary row, showing '--' for unavailable values
%   (published tables do not quote every metric).
s = sprintf('%-20s', label);
for k = 1:numel(v)
    if isnan(v(k))
        s = [s sprintf('%10s','--')]; %#ok<AGROW>
    else
        s = [s sprintf('%10.1f', v(k))]; %#ok<AGROW>
    end
end
fprintf('%s\n', s);
end

function localPrintRestSummary(rest)
a = localFirstFinite(rest.APD90);
b = localLastFinite(rest.APD90);
fprintf('%-28s %10.1f %10.1f   (%d S1 beats, %d DIs tested)\n', ...
    rest.modelName, a, b, rest.nS1, numel(rest.DI));
end

function x = localFirstFinite(v)
idx = find(isfinite(v), 1, 'first');
if isempty(idx), x = NaN; else, x = v(idx); end
end

function x = localLastFinite(v)
idx = find(isfinite(v), 1, 'last');
if isempty(idx), x = NaN; else, x = v(idx); end
end

function row = localMetricRow(model, config, m, ref, src)
row = {model, config, m.RMP, m.overshoot, m.amplitude, m.plateau, m.dVdtMax, ...
       m.APD25, m.APD50, m.APD90, m.CaiRest, m.CaiPeak, ...
       ref.RMP, ref.plateau, ref.dVdtMax, ref.APD25, ref.APD50, ref.APD90, ...
       ref.CaiRest, ref.CaiPeak, src};
end

function [tRel, idx] = localLastBeatWindow(res, ap, durMs)
%LOCALLASTBEATWINDOW Index range covering the last paced beat, with time
%   re-zeroed to that beat's stimulus so beats from different runs (which
%   may have different APDs and even different time steps) can be overlaid.
tStim = ap.stimTimes(ap.count);
idx   = find(res.t >= tStim & res.t <= tStim + durMs);
tRel  = res.t(idx) - tStim;
end

function fig = localPlotFailingOverlay(resNF, resF, apNF, apF, pal)
%LOCALPLOTFAILINGOVERLAY Failing vs nonfailing action potential and
%   [Ca2+]i transient -- the central result of Priebe & Beuckelmann (1998).
%   Colours: the nonfailing trace takes the palette's neutral V colour and
%   the failing trace the INaCa colour, because enhanced Na+/Ca2+ exchange
%   is the main driver of the prolongation being shown.
DUR = 700;   % ms of the beat to show
[tN, iN] = localLastBeatWindow(resNF, apNF, DUR);
[tF, iF] = localLastBeatWindow(resF,  apF,  DUR);
mN = localLastBeatMetrics(resNF, apNF);
mF = localLastBeatMetrics(resF,  apF);

fig = figure('Color','w','Position',[100 100 900 780]);
tl = tiledlayout(2,1,'TileSpacing','compact');

nexttile;
plot(tN, resNF.V(iN), 'Color', pal.V,     'LineWidth', 1.8, 'DisplayName', 'Nonfailing'); hold on;
plot(tF, resF.V(iF),  'Color', pal.INaCa, 'LineWidth', 1.8, 'DisplayName', 'Failing');
% Properties set after construction rather than as trailing arguments:
% YLINE would otherwise be free to read a trailing char as a line label.
hZero = yline(0, ':');
hZero.Color = pal.grid;  hZero.HandleVisibility = 'off';
grid on; box off;
ylabel('V_m (mV)');
title('Action potential prolongation in heart failure');
legend('Location','northeast');
text(0.42, 0.62, sprintf(['APD_{90}\n  nonfailing %.0f ms (paper 374)\n' ...
    '  failing    %.0f ms (paper 549)\n  prolongation %+.0f%%'], ...
    mN.APD90, mF.APD90, 100*(mF.APD90/mN.APD90-1)), ...
    'Units','normalized', 'FontSize', 9, 'VerticalAlignment','top');

nexttile;
plot(tN, resNF.Cai(iN)*1e6, 'Color', pal.V,     'LineWidth', 1.8, 'DisplayName', 'Nonfailing'); hold on;
plot(tF, resF.Cai(iF)*1e6,  'Color', pal.INaCa, 'LineWidth', 1.8, 'DisplayName', 'Failing');
grid on; box off;
xlabel('Time after stimulus (ms)');
ylabel('[Ca^{2+}]_i (nM)');
title('Ca^{2+} transient: smaller peak and slower decay in failure');
legend('Location','northeast');

title(tl, {'Priebe & Beuckelmann (1998): failing vs. nonfailing human ventricular myocyte', ...
    'last beat of a 1 Hz train; prolongation is driven by enhanced I_{NaCa}, slowed [Ca^{2+}]_i decay, reduced I_{K1} and I_{NaK}'});
end

function fig = localPlotCellTypes(resList, mList, apList, labels, pal)
%LOCALPLOTCELLTYPES Transmural heterogeneity in the reduced model --
%   epicardial, midmyocardial and endocardial action potentials, mirroring
%   Bernus et al. (2002) Fig. 6A-C. The palette has no cell-type entries,
%   so three well-separated existing palette colours are reused rather
%   than hardcoding new RGB triples (house rule).
cols = {pal.reduced, pal.f, pal.IK1};
DUR = 600;

fig = figure('Color','w','Position',[100 100 1050 620]);
tl = tiledlayout(1,2,'TileSpacing','compact');

nexttile;
hold on;
for k = 1:numel(resList)
    [tRel, idx] = localLastBeatWindow(resList{k}, apList{k}, DUR);
    plot(tRel, resList{k}.V(idx), 'Color', cols{k}, 'LineWidth', 1.8, 'DisplayName', labels{k});
end
grid on; box off;
xlabel('Time after stimulus (ms)'); ylabel('V_m (mV)');
title('Action potentials by cell type');
legend('Location','northeast');

nexttile;
apd50 = cellfun(@(m) m.APD50, mList);
apd90 = cellfun(@(m) m.APD90, mList);
b = bar(1:numel(labels), [apd50(:) apd90(:)]);
b(1).FaceColor = pal.to;   b(1).DisplayName = 'APD_{50}';
b(2).FaceColor = pal.X;    b(2).DisplayName = 'APD_{90}';
xticks(1:numel(labels)); xticklabels({'Epi','M','Endo'});
grid on; box off;
ylabel('Duration (ms)');
title('APD by cell type');
legend('Location','northwest');

title(tl, {'Bernus et al. (2002) reduced model: transmural heterogeneity (cf. Fig. 6A-C)', ...
    'epi / M / endo differ only in g_{to}, g_K and the I_{to} kinetics (1/p, V_{shift})'});
end

function fig = localPlotRestitution(restFull, restRed, pal)
%LOCALPLOTRESTITUTION S1-S2 APD90 restitution curves for both models.
%   Points where the S2 failed to capture (or failed to repolarize inside
%   the record) are simply absent -- COMPUTE_PB_RESTITUTION returns NaN
%   there, which MATLAB skips, so a gap in a curve means "no capture".
fig = figure('Color','w','Position',[100 100 900 640]);
plot(restFull.DI, restFull.APD90, '-o', 'Color', pal.full, ...
    'MarkerFaceColor', pal.full, 'LineWidth', 1.8, ...
    'DisplayName', sprintf('%s, %d S1 beats', restFull.modelName, restFull.nS1)); hold on;
plot(restRed.DI, restRed.APD90, '-s', 'Color', pal.reduced, ...
    'MarkerFaceColor', pal.reduced, 'LineWidth', 1.8, ...
    'DisplayName', sprintf('%s, %d S1 beats', restRed.modelName, restRed.nS1));
grid on; box off;
xlabel('Diastolic interval (ms)');
ylabel('S2 APD_{90} (ms)');
title({'APD restitution, S1-S2 protocol (Bernus et al. 2002)', ...
    'Bernus reports the reduced model matching the Morgan et al. (1992) human RV data slightly better than the full model'});
legend('Location','southeast');
text(0.02, 0.05, ...
    sprintf(['Full-model curve uses only %d S1 beats for runtime reasons;\n' ...
             'per Bernus the full model needs >1 min of model time to\n' ...
             're-equilibrate, so it is further from steady state here.'], restFull.nS1), ...
    'Units','normalized', 'FontSize', 8, 'VerticalAlignment','bottom');
end
