function s = detect_pb_ap(res, opts)
arguments
    res (1,1) struct
    opts.repolFractions (1,:) double = [0.25 0.5 0.90]
    opts.minISI (1,1) double = 200
    opts.searchWindow (1,1) double = 900
    % Capture criteria: a beat must show a genuine regenerative upstroke
    % before any APD is reported for it. Defaults are deliberately loose
    opts.minAmplitude (1,1) double = 60
    opts.minOvershoot (1,1) double = 0
end

t = res.t; V = res.V; Istim = res.Istim;
dt = t(2)-t(1);

stimOnIdx = find(Istim(1:end-1) <= 0 & Istim(2:end) > 0);
if ~isempty(Istim) && Istim(1) > 0
    stimOnIdx = [1, stimOnIdx];
end
keep = true(size(stimOnIdx));
last = 1;
for k = 2:numel(stimOnIdx)
    if t(stimOnIdx(k)) - t(stimOnIdx(last)) < opts.minISI
        keep(k) = false;
    else
        last = k;
    end
end
stimOnIdx = stimOnIdx(keep);
nBeats = numel(stimOnIdx);

if nBeats == 0
    s = struct('stimTimes',[], 'count',0, 'RMP',[], 'overshoot',[], 'amplitude',[], ...
               'dVdtMax',[], 'APD',zeros(0,numel(opts.repolFractions)), ...
               'captured',true(1,0), 'repolFractions',opts.repolFractions);
    return
end

winSamp = round(opts.searchWindow/dt);
RMP = zeros(1,nBeats); overshoot = zeros(1,nBeats); dVdtMax = zeros(1,nBeats);
APD = nan(nBeats, numel(opts.repolFractions));
captured = false(1,nBeats);
dVdt_full = [0, diff(V)./diff(t)];

for k = 1:nBeats
    iStart = stimOnIdx(k);
    if k < nBeats
        iEnd = min([numel(V), iStart+winSamp, stimOnIdx(k+1)-1]);
    else
        iEnd = min(numel(V), iStart+winSamp);
    end
    RMP(k) = V(iStart);
    [pk, relIdx] = max(V(iStart:iEnd));
    overshoot(k) = pk;
    pkIdx = iStart+relIdx-1;
    dVdtMax(k) = max(dVdt_full(iStart:pkIdx));
    amp = overshoot(k) - RMP(k);

    captured(k) = (amp >= opts.minAmplitude) && (overshoot(k) >= opts.minOvershoot);
    if ~captured(k)
        RMP(k) = V(iStart);          % keep: still a valid resting measurement
        overshoot(k) = pk;           % keep: reports how far it actually got
        dVdtMax(k) = NaN;
        continue                     % leave APD(k,:) as NaN
    end

    seg = V(pkIdx:iEnd); segT = t(pkIdx:iEnd);
    for j = 1:numel(opts.repolFractions)
        target = overshoot(k) - opts.repolFractions(j)*amp;
        fallIdx = find(seg(1:end-1) >= target & seg(2:end) < target, 1, 'last');
        if ~isempty(fallIdx)
            f = (seg(fallIdx)-target)/(seg(fallIdx)-seg(fallIdx+1));
            tRepol = segT(fallIdx) + f*(segT(fallIdx+1)-segT(fallIdx));
            APD(k,j) = tRepol - t(iStart);
        end
    end
end

s.stimTimes = t(stimOnIdx);
s.count = nBeats;
s.RMP = RMP;
s.overshoot = overshoot;
s.amplitude = overshoot - RMP;
s.dVdtMax = dVdtMax;
s.APD = APD;
s.captured = captured;   
s.repolFractions = opts.repolFractions;
end
