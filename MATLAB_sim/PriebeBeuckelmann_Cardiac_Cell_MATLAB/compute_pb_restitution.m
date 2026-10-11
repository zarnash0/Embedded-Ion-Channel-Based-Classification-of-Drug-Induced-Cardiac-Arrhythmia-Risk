function rest = compute_pb_restitution(modelFcn, params, opts)
%   Name-value options:
%     nS1        number of S1 conditioning stimuli                  [10]
%     S1period   S1 basic cycle length (ms)                       [1000]
%     DIlist     diastolic intervals to test (ms)
%                                 [20 30 50 75 100 150 200 300 400 500]
%     dt         integration step (ms); [] auto-selects 0.005 for the
%                full model and 0.02 for the reduced one             []
%     pulseDur   stimulus pulse width (ms)                           [3]
%     amplitude  stimulus amplitude (pA/pF)                         [15]
%     startTime  time of the first S1 (ms)                          [50]
%     postS2     record length kept after the S2 (ms)              [800]
%     searchWindow, minISI   passed through to DETECT_PB_AP  [1000, 100]
%     captureThreshold  peak V (mV) below which the S2 is treated as
%                having failed to capture (cell still refractory)    [0]
%     verbose    print per-DI progress                            [true]

arguments
    modelFcn (1,1) function_handle
    params (1,1) struct
    opts.nS1 (1,1) double {mustBePositive, mustBeInteger} = 10
    opts.S1period (1,1) double {mustBePositive} = 1000
    opts.DIlist (1,:) double {mustBePositive} = [20 30 50 75 100 150 200 300 400 500]
    opts.dt (1,:) double {mustBePositive} = []
    opts.pulseDur (1,1) double {mustBePositive} = 3
    opts.amplitude (1,1) double = 15
    opts.startTime (1,1) double {mustBeNonnegative} = 50
    opts.postS2 (1,1) double {mustBePositive} = 800
    opts.searchWindow (1,1) double {mustBePositive} = 1000
    opts.minISI (1,1) double {mustBePositive} = 100
    opts.captureThreshold (1,1) double = 0
    opts.verbose (1,1) logical = true
end

modelStr = func2str(modelFcn);
if isempty(opts.dt)
    if contains(lower(modelStr), 'full')
        dt = 0.005;
    elseif contains(lower(modelStr), 'reduced')
        dt = 0.02;
    else
        dt = 0.005;   
        warning('compute_pb_restitution:dtGuess', ...
            ['Could not identify the model from the handle "%s"; using the ' ...
             'conservative dt = %.4g ms. Pass dt explicitly to override.'], modelStr, dt);
    end
else
    dt = opts.dt(1);
end

S1times  = opts.startTime + (0:opts.nS1-1)*opts.S1period;
tLastS1  = S1times(end);

%% 
tEndS1  = tLastS1 + opts.S1period;
tGridS1 = 0:dt:tEndS1;
stimS1  = pb_rest_stim_vector(tGridS1, S1times, opts.pulseDur, opts.amplitude);

resS1 = modelFcn(params, tEndS1, 'dt', dt, 'stim', stimS1);
apS1  = detect_pb_ap(resS1, 'repolFractions', 0.90, ...
                     'minISI', opts.minISI, 'searchWindow', opts.searchWindow);

if apS1.count == 0
    error('compute_pb_restitution:noS1beats', ...
        'No S1 beats were detected in the reference run -- check the stimulus settings.');
end
lastGood = find(~isnan(apS1.APD(:,1)), 1, 'last');
if isempty(lastGood)
    error('compute_pb_restitution:noS1repol', ...
        ['No S1 beat repolarized to 90%% within the %.0f ms measurement window; ' ...
         'increase searchWindow or S1period.'], opts.searchWindow);
end
APD90_S1 = apS1.APD(lastGood,1);

modelName = resS1.model;
if isfield(params, 'condition')
    modelName = sprintf('%s (%s)', modelName, params.condition);
elseif isfield(params, 'cellType')
    modelName = sprintf('%s (%s)', modelName, params.cellType);
end

if opts.verbose
    fprintf('  %s restitution: %d S1 at %.0f ms, dt = %.4g ms, last S1 APD90 = %.1f ms\n', ...
        modelName, opts.nS1, opts.S1period, dt, APD90_S1);
end

%% 
nDI          = numel(opts.DIlist);
APD90        = nan(1, nDI);
S2times      = nan(1, nDI);
S2overshoot  = nan(1, nDI);
captured     = false(1, nDI);

for k = 1:nDI
    DI  = opts.DIlist(k);
    tS2 = tLastS1 + APD90_S1 + DI;
    tEnd   = tS2 + opts.postS2;
    tGrid  = 0:dt:tEnd;
    stim   = pb_rest_stim_vector(tGrid, [S1times tS2], opts.pulseDur, opts.amplitude);

    res = modelFcn(params, tEnd, 'dt', dt, 'stim', stim);
    ap  = detect_pb_ap(res, 'repolFractions', 0.90, ...
                       'minISI', opts.minISI, 'searchWindow', opts.searchWindow);
    if ap.count == 0 || abs(ap.stimTimes(end) - tS2) > max(2*opts.pulseDur, 2*dt)
        warning('compute_pb_restitution:s2NotDetected', ...
            'S2 at t = %.1f ms (DI = %.1f ms) was not resolved as a separate beat; skipping.', tS2, DI);
        continue
    end

    S2times(k)     = ap.stimTimes(end);
    S2overshoot(k) = ap.overshoot(end);
    captured(k)    = ap.captured(end) && (ap.overshoot(end) >= opts.captureThreshold);
    if captured(k)
        APD90(k) = ap.APD(end,1);   
    end

    if opts.verbose
        if ~captured(k)
            fprintf('    DI = %6.1f ms -> no capture (peak %.1f mV, cell still refractory)\n', ...
                DI, S2overshoot(k));
        elseif isnan(APD90(k))
            fprintf('    DI = %6.1f ms -> captured but did not repolarize within the window\n', DI);
        else
            fprintf('    DI = %6.1f ms -> S2 APD90 = %6.1f ms\n', DI, APD90(k));
        end
    end
end

rest = struct( ...
    'DI',          opts.DIlist, ...
    'APD90',       APD90, ...
    'modelName',   modelName, ...
    'S1_APD90',    APD90_S1, ...
    'S2times',     S2times, ...
    'S2overshoot', S2overshoot, ...
    'captured',    captured, ...
    'nS1',         opts.nS1, ...
    'S1period',    opts.S1period, ...
    'dt',          dt);
end

% 
function I = pb_rest_stim_vector(t, onsets, pulseDur, amplitude)
I = zeros(size(t));
for k = 1:numel(onsets)
    I(t >= onsets(k) & t < onsets(k) + pulseDur) = amplitude;
end
end
