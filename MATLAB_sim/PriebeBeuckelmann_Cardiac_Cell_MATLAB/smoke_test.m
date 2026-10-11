%% SMOKE_TEST - exercise every plotting/analysis function on tiny runs.
% Purpose: catch runtime errors (bad field names, wrong argument counts,
% empty-beat edge cases) in seconds, instead of discovering them 40
% minutes into RUN_PB_DEMO. Figures are created invisibly and closed.
clear; clc; close all;
% onCleanup so an error mid-test cannot leave the session with figures
% permanently invisible.
set(0,'DefaultFigureVisible','off');
restoreFigs = onCleanup(@() set(0,'DefaultFigureVisible','on'));
ok = @(s) fprintf('  ok  %s\n', s);
fail = 0;

T = 1500;   % just long enough for one full beat plus baseline

fprintf('Simulating (short)...\n');
pFull = pb_full_params("nonfailing");
pRed  = pb_reduced_params("epi");
rFull = simulate_PB_full(pFull, T, dt=0.005);
rRed  = simulate_PB_reduced(pRed, T, dt=0.02);
fprintf('  full: %d beats | reduced: %d beats\n', ...
    detect_pb_ap(rFull).count, detect_pb_ap(rRed).count);

tests = {
  'plot_pb_ap(full)',        @() plot_pb_ap(rFull)
  'plot_pb_ap(reduced)',     @() plot_pb_ap(rRed)
  'plot_pb_currents(full)',  @() plot_pb_currents(rFull)
  'plot_pb_currents(red)',   @() plot_pb_currents(rRed)
  'plot_pb_gates(full)',     @() plot_pb_gates(rFull)
  'plot_pb_gates(reduced)',  @() plot_pb_gates(rRed)
  'plot_pb_calcium(full)',   @() plot_pb_calcium(rFull)
  'plot_pb_phase3d(full)',   @() plot_pb_phase3d(rFull)
  'plot_pb_phase3d(red)',    @() plot_pb_phase3d(rRed)
};
for k = 1:size(tests,1)
    try
        f = tests{k,2}(); ok(tests{k,1}); close(f);
    catch ME
        fprintf('  FAIL %s\n       %s\n', tests{k,1}, ME.message);
        fail = fail + 1;
    end
end

% plot_pb_calcium must REFUSE the reduced model, with a clear message.
try
    plot_pb_calcium(rRed);
    fprintf('  FAIL plot_pb_calcium(reduced) should have errored\n'); fail = fail+1;
catch
    ok('plot_pb_calcium(reduced) correctly rejected');
end

fprintf('Comparison (short)...\n');
try
    cmp = compare_pb_models(tEnd=T, dtFull=0.005, dtReduced=0.02, verbose=false);
    ok('compare_pb_models');
    f = plot_pb_comparison(cmp); ok('plot_pb_comparison'); close(f);
    fprintf('    speedup=%.2fx  sameDt=%.2fx  states %d vs %d\n', ...
        cmp.speedup, cmp.speedupSameDt, cmp.nStateFull, cmp.nStateReduced);
catch ME
    fprintf('  FAIL comparison\n       %s\n', ME.message);
    fprintf('       at %s line %d\n', ME.stack(1).name, ME.stack(1).line);
    fail = fail + 1;
end

fprintf('Restitution (2 DIs, reduced only)...\n');
try
    rr = compute_pb_restitution(@simulate_PB_reduced, pRed, ...
        nS1=3, DIlist=[100 400], dt=0.02, verbose=false);
    ok(sprintf('compute_pb_restitution -> APD90 = [%s]', num2str(rr.APD90,'%.1f ')));
catch ME
    fprintf('  FAIL restitution\n       %s\n', ME.message);
    fprintf('       at %s line %d\n', ME.stack(1).name, ME.stack(1).line);
    fail = fail + 1;
end

% --- Regression guard: the optimized RHS must match the reference ------
fprintf('Fast-path equivalence (both parameter conditions)...\n');
try
    verify_pb_rhs(true);          % raises on mismatch
    ok('verify_pb_rhs');
catch ME
    fprintf('  FAIL verify_pb_rhs\n       %s\n', ME.message); fail = fail + 1;
end

% --- Configurations not covered above ----------------------------------
fprintf('Other configurations (failing / M / endo / rk4)...\n');
cfgs = {
  'full failing',  @() simulate_PB_full(pb_full_params("failing"), T, dt=0.005)
  'reduced M',     @() simulate_PB_reduced(pb_reduced_params("M"),    T, dt=0.02)
  'reduced endo',  @() simulate_PB_reduced(pb_reduced_params("endo"), T, dt=0.02)
  'reduced rk4',   @() simulate_PB_reduced(pb_reduced_params("epi"),  T, dt=0.005, method="rk4")
};
for k = 1:size(cfgs,1)
    try
        r = cfgs{k,2}();
        a = detect_pb_ap(r);
        kk = find(all(isfinite(a.APD),2)' & a.captured, 1, 'last');
        % Coarse physiological assertions: these would catch a gross
        % regression (a dead cell, a runaway, a sign flip) without being
        % so tight that a legitimate parameter change trips them.
        assert(~isempty(kk), 'no measurable beat');
        assert(a.RMP(kk) > -100 && a.RMP(kk) < -80, 'RMP %.1f out of range', a.RMP(kk));
        assert(a.APD(kk,3) > 150 && a.APD(kk,3) < 800, 'APD90 %.1f out of range', a.APD(kk,3));
        assert(a.dVdtMax(kk) > 100, 'dV/dt max %.1f too low', a.dVdtMax(kk));
        ok(sprintf('%-13s RMP=%.1f APD90=%.1f dVdt=%.0f', cfgs{k,1}, ...
            a.RMP(kk), a.APD(kk,3), a.dVdtMax(kk)));
    catch ME
        fprintf('  FAIL %s\n       %s\n', cfgs{k,1}, ME.message); fail = fail + 1;
    end
end

% --- Non-capture must be reported as such, not as a tiny APD -----------
try
    weak = simulate_PB_reduced(pb_reduced_params("epi"), 1200, dt=0.02, ...
        stim=pb_make_train(1000, 3, 0.5, 50));    % far sub-threshold
    aw = detect_pb_ap(weak);
    assert(~any(aw.captured), 'sub-threshold stimulus was reported as captured');
    assert(all(isnan(aw.APD(:))), 'sub-threshold beat returned a finite APD');
    ok('non-capture correctly flagged (captured=false, APD=NaN)');
catch ME
    fprintf('  FAIL non-capture guard\n       %s\n', ME.message); fail = fail + 1;
end

fprintf('\n%s (%d failures)\n', repmat('=',1,40), fail);
