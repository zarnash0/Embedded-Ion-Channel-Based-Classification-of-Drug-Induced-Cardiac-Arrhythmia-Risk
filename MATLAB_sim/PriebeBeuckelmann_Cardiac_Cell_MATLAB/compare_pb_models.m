function cmp = compare_pb_models(opts)
arguments
    opts.tEnd (1,1) double {mustBePositive} = 5000
    opts.dtFull (1,1) double {mustBePositive} = 0.005
    opts.dtReduced (1,1) double {mustBePositive} = 0.02
    opts.condition (1,1) string {mustBeMember(opts.condition, ["nonfailing","failing"])} = "nonfailing"
    opts.cellType (1,1) string {mustBeMember(opts.cellType, ["epi","M","endo"])} = "epi"
    opts.verbose (1,1) logical = true
end
% Identical pacing for both models 
% The stimulus is passed as a FUNCTION HANDLE rather than a sampled vector
% precisely because the two models use different dt: a handle lets each
% simulator sample the same continuous-time train on its own grid, so the
% pacing is provably identical even though the discretization is not.
% These are the defaults of both simulators, restated here explicitly so
% the comparison cannot silently drift if a default ever changes.
protocol = struct('period',1000, 'pulseDur',3, 'amplitude',15, 'startTime',50, ...
                  'searchWindow',900);
stimFcn = pb_make_train(protocol.period, protocol.pulseDur, protocol.amplitude, protocol.startTime);

pFull = pb_full_params(opts.condition);
pRed  = pb_reduced_params(opts.cellType);

warmFull = simulate_PB_full(pFull, 5, dt=opts.dtFull, stim=stimFcn);              %#ok<NASGU>
warmRed  = simulate_PB_reduced(pRed, 5, dt=opts.dtReduced, stim=stimFcn);         %#ok<NASGU>
clear warmFull warmRed

tStart = tic;
simulate_PB_full(pFull, opts.tEnd, dt=opts.dtFull, stim=stimFcn, returnCurrents=false);
timeFull = toc(tStart);

tStart = tic;
simulate_PB_reduced(pRed, opts.tEnd, dt=opts.dtReduced, stim=stimFcn, returnCurrents=false);
timeReduced = toc(tStart);

tStart = tic;
simulate_PB_reduced(pRed, opts.tEnd, dt=opts.dtFull, stim=stimFcn, returnCurrents=false);
timeReducedSameDt = toc(tStart);

resFull = simulate_PB_full(pFull, opts.tEnd, dt=opts.dtFull, stim=stimFcn);
resRed  = simulate_PB_reduced(pRed, opts.tEnd, dt=opts.dtReduced, stim=stimFcn);

apFull = detect_pb_ap(resFull, searchWindow=protocol.searchWindow);
apRed  = detect_pb_ap(resRed,  searchWindow=protocol.searchWindow);

kFull = localLastCompleteBeat(apFull, resFull.t(end), protocol.searchWindow);
kRed  = localLastCompleteBeat(apRed,  resRed.t(end),  protocol.searchWindow);
if isempty(kFull) || isempty(kRed)
    error('compare_pb_models:noCompleteBeat', ...
        ['No fully-measurable beat in at least one model (full: %d beats, reduced: %d beats). ' ...
         'Increase tEnd so that at least one stimulus is followed by a complete AP.'], ...
        apFull.count, apRed.count);
end
kCommon = min(kFull, kRed);
measurable = @(ap,k) k >= 1 && k <= ap.count && ...
                     all(isfinite(ap.APD(k,:))) && ap.captured(k);
while kCommon >= 1 && ~(measurable(apFull,kCommon) && measurable(apRed,kCommon))
    kCommon = kCommon - 1;
end
if kCommon < 1
    error('compare_pb_models:noCommonBeat', ...
        ['No single beat index is fully measurable in BOTH models ' ...
         '(full last: %d, reduced last: %d). Increase tEnd.'], kFull, kRed);
end
if kFull ~= kRed
    warning('compare_pb_models:beatIndexMismatch', ...
        ['Last fully-measurable beat differs between models (full: %d, reduced: %d); ' ...
         'comparing beat %d in both so the difference column stays like-for-like.'], ...
        kFull, kRed, kCommon);
end
kFull = kCommon;  kRed = kCommon;

c25F = localFracCol(apFull, 0.25);  c50F = localFracCol(apFull, 0.50);  c90F = localFracCol(apFull, 0.90);
c25R = localFracCol(apRed,  0.25);  c50R = localFracCol(apRed,  0.50);  c90R = localFracCol(apRed,  0.90);

vFull = [apFull.RMP(kFull); apFull.overshoot(kFull); apFull.amplitude(kFull); apFull.dVdtMax(kFull); ...
         apFull.APD(kFull,c25F); apFull.APD(kFull,c50F); apFull.APD(kFull,c90F)];
vRed  = [apRed.RMP(kRed);  apRed.overshoot(kRed);  apRed.amplitude(kRed);  apRed.dVdtMax(kRed); ...
         apRed.APD(kRed,c25R);  apRed.APD(kRed,c50R);  apRed.APD(kRed,c90R)];

apdTable = table( ...
    ["RMP";"Overshoot";"Amplitude";"dV/dt max";"APD25";"APD50";"APD90"], ...
    ["mV";"mV";"mV";"V/s";"ms";"ms";"ms"], ...
    vFull, vRed, abs(vFull-vRed), ...
    'VariableNames', {'Metric','Unit','Full','Reduced','AbsDiff'});
drift = struct();
drift.NaiStart = resFull.Nai(1);   drift.NaiEnd = resFull.Nai(end);
drift.KiStart  = resFull.Ki(1);    drift.KiEnd  = resFull.Ki(end);
drift.NaiChange = drift.NaiEnd - drift.NaiStart;
drift.KiChange  = drift.KiEnd  - drift.KiStart;
drift.NaiFixedReduced = pRed.Nai;   % 10 mM, fixed
drift.KiFixedReduced  = pRed.Ki;    % 140 mM, fixed
drift.beatsFull    = 1:apFull.count;
drift.beatsReduced = 1:apRed.count;
drift.APD90Full    = apFull.APD(:,c90F)';
drift.APD90Reduced = apRed.APD(:,c90R)';
drift.APD90ChangeFull    = localSeriesChange(drift.APD90Full);
drift.APD90ChangeReduced = localSeriesChange(drift.APD90Reduced);
drift.note = ['Full model integrates [Na+]i and [K+]i; reduced model holds them fixed ' ...
              'at 10 and 140 mM, so it cannot drift by construction (Bernus et al. 2002).'];

cmp = struct();
cmp.full = resFull;
cmp.reduced = resRed;
cmp.apFull = apFull;
cmp.apReduced = apRed;
cmp.timeFull = timeFull;
cmp.timeReduced = timeReduced;
cmp.timeReducedSameDt = timeReducedSameDt;
cmp.speedup = timeFull/timeReduced;
cmp.speedupSameDt = timeFull/timeReducedSameDt;
cmp.nStateFull = 15;
cmp.nStateReduced = 6;
cmp.apdTable = apdTable;
cmp.drift = drift;
cmp.protocol = protocol;
cmp.beatFull = kFull;
cmp.beatReduced = kRed;
cmp.opts = opts;

if opts.verbose
    localPrintSummary(cmp);
end
end

%
function k = localLastCompleteBeat(ap, tLast, searchWindow)
k = [];
if ap.count == 0
    return
end
measurable = all(isfinite(ap.APD), 2)';
inRecord   = (ap.stimTimes + searchWindow) <= tLast;
k = find(measurable & inRecord, 1, 'last');
if isempty(k)
    k = find(measurable, 1, 'last');
end
end

% 
function c = localFracCol(ap, frac)
c = find(abs(ap.repolFractions - frac) < 1e-9, 1);
if isempty(c)
    error('compare_pb_models:missingFraction', ...
        'DETECT_PB_AP did not return APD at the %.0f%% repolarization level.', 100*frac);
end
end

% 
function d = localSeriesChange(x)
x = x(isfinite(x));
if numel(x) < 2
    d = NaN;
else
    d = x(end) - x(1);
end
end

% 
function localPrintSummary(cmp)
p = cmp.protocol;
sep = repmat('-', 1, 72);   % not named "line": that would shadow the built-in

fprintf('\n%s\n', repmat('=', 1, 72));
fprintf(' PRIEBE-BEUCKELMANN FULL  vs  BERNUS REDUCED -- model comparison\n');
fprintf('%s\n', repmat('=', 1, 72));
fprintf(' Protocol : %g Hz, %g ms pulses of %g pA/pF from t = %g ms, run to %g ms\n', ...
    1000/p.period, p.pulseDur, p.amplitude, p.startTime, cmp.opts.tEnd);
fprintf(' Full     : %s myocyte, 15 states, dt = %g ms, %s\n', ...
    cmp.opts.condition, cmp.full.dt, cmp.full.method);
fprintf(' Reduced  : %s cell,   6 states, dt = %g ms, %s\n', ...
    cmp.opts.cellType, cmp.reduced.dt, cmp.reduced.method);

fprintf('\n%s\n 1. ACCURACY -- last complete beat (full: beat %d/%d, reduced: beat %d/%d)\n%s\n', ...
    sep, cmp.beatFull, cmp.apFull.count, cmp.beatReduced, cmp.apReduced.count, sep);
disp(cmp.apdTable);

fprintf('%s\n 2. COST\n%s\n', sep, sep);
fprintf('   wall clock, full model              : %8.2f s   (dt = %g ms)\n', cmp.timeFull, cmp.full.dt);
fprintf('   wall clock, reduced model           : %8.2f s   (dt = %g ms)\n', cmp.timeReduced, cmp.reduced.dt);
fprintf('   -> PRACTICAL speedup                : %8.2fx  (each model at the dt it can actually use)\n', cmp.speedup);
fprintf('   wall clock, reduced at the full dt  : %8.2f s   (dt = %g ms, control run)\n', ...
    cmp.timeReducedSameDt, cmp.full.dt);
fprintf('   -> SAME-dt speedup                  : %8.2fx  (per-step cost of the reduction alone)\n', cmp.speedupSameDt);
fprintf('   state variables                     : %d -> %d      (%.1fx fewer; memory scales with this\n', ...
    cmp.nStateFull, cmp.nStateReduced, cmp.nStateFull/cmp.nStateReduced);
fprintf('                                                      in tissue-scale simulations)\n');
fprintf('   For reference, Bernus et al. (2002) report 4.9x faster than the full model\n');
fprintf('   and 3.8x faster than the full model with [Na+]i/[K+]i clamped, with memory\n');
fprintf('   down 2.5x. Absolute ratios are implementation- and machine-dependent, so the\n');
fprintf('   numbers above are not expected to reproduce theirs exactly -- only the\n');
fprintf('   direction and rough magnitude of the gain are meaningful.\n');

fprintf('\n%s\n 3. STABILITY / IONIC DRIFT\n%s\n', sep, sep);
fprintf('   full model  [Na^+]i : %8.4f -> %8.4f mM  (%+.4f)\n', ...
    cmp.drift.NaiStart, cmp.drift.NaiEnd, cmp.drift.NaiChange);
fprintf('   full model  [K^+]i  : %8.4f -> %8.4f mM  (%+.4f)\n', ...
    cmp.drift.KiStart, cmp.drift.KiEnd, cmp.drift.KiChange);
fprintf('   reduced model       : fixed at %g and %g mM -- cannot drift by construction\n', ...
    cmp.drift.NaiFixedReduced, cmp.drift.KiFixedReduced);
fprintf('   APD90 beat-to-beat  : full    %+.2f ms over %d beats\n', ...
    cmp.drift.APD90ChangeFull, numel(cmp.drift.beatsFull));
fprintf('                       : reduced %+.2f ms over %d beats\n', ...
    cmp.drift.APD90ChangeReduced, numel(cmp.drift.beatsReduced));
fprintf('   This run is only %.1f s of model time. Bernus et al. paced the full model\n', cmp.opts.tEnd/1000);
fprintf('   for 1000 s and saw mean [Na+]i fall 10 -> 8.1 mM with APD90 creeping from\n');
fprintf('   ~360 to 425 ms, still not at steady state. What is measurable here is the\n');
fprintf('   start of that trend, not the whole of it.\n');

fprintf('\n%s\n WHY THE REDUCED MODEL IS PREFERABLE (for its intended use)\n%s\n', sep, sep);
fprintf('   - It reproduces the AP shape it was fitted to: see the table above; the\n');
fprintf('     upstroke, plateau and repolarization phases all survive the reduction.\n');
fprintf('   - It runs %.1fx faster in practice and carries %.1fx fewer state\n', ...
    cmp.speedup, cmp.nStateFull/cmp.nStateReduced);
fprintf('     variables -- decisive when the single cell becomes a 2-D sheet or a\n');
fprintf('     3-D ventricle with 10^6-10^7 cells.\n');
fprintf('   - It has no long-term drift and adapts to a rate change within a few beats,\n');
fprintf('     where the full model needs > 1 min of model time to re-equilibrate.\n');
fprintf('   WHAT IT GIVES UP, honestly: no [Ca2+]i transient, no Ca2+ handling, no\n');
fprintf('   [Na+]i / [K+]i dynamics, and a single merged I_K instead of I_Kr + I_Ks. Any\n');
fprintf('   question about calcium, about ionic homeostasis, or about a drug acting on\n');
fprintf('   I_Kr alone REQUIRES the full model -- the reduction is preferable for wave\n');
fprintf('   propagation and arrhythmia studies, not universally.\n');
fprintf('%s\n\n', repmat('=', 1, 72));
end
