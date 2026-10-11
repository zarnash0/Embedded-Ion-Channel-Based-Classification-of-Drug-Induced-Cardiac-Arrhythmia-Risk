function s = detect_beats(res, opts)
%DETECT_BEATS Beat detector and electrophysiological analytics for a Noble-model run.
%   takes the struct returned by SIMULATE_NOBLE and
%   finds each spontaneous upstroke then measures, per beat:
%     overshoot        peak (mV) reached during the upstroke
%     MDP              maximum diastolic potential (mV) -- the minimum V
%                      in the diastolic interval immediately preceding
%                      this beat's upstroke
%     amplitude        overshoot - MDP (mV)
%     dVdtMax          maximum upstroke velocity (mV/ms == V/s, the
%                      conventional cardiac electrophysiology unit)
%     APD              action-potential duration (ms) at
%                      opts.repolFraction repolarization (default APD90)

%   Also computes two quantitative checks against the specific,
%   citable claims in Fink & Noble:
%     pacemakerMechanismScore    percentage of samples in the diastolic
%       interval preceding each beat where n is non-increasing AND V is
%       non-decreasing simultaneously -- i.e. how consistent the
%       simulated diastolic depolarization is with this model's actual
%       mechanism (decaying i_K2 permitting net inward current), as
%       opposed to some other/incorrect dynamic. Should be high (>90%).
%     persistentNaConductancePct    ratio (%) of the *effective* total
%       Na+ conductance (gNaMax*m^3*h + gNaBg) during the plateau to its
%       peak value at the upstroke -- Scholarpedia's Fig. 4 caption
%       states this should be "about 2% of the peak conductance."

arguments
    res (1,1) struct
    opts.threshold (1,1) double = -40
    opts.refractory (1,1) double = 150
    opts.peakWindow (1,1) double = 20
    opts.repolFraction (1,1) double = 0.9
end

t = res.t; V = res.V; n = res.n; gNa = res.gNa;
dt = t(2) - t(1);

crossIdx = find(V(1:end-1) < opts.threshold & V(2:end) >= opts.threshold);
if isempty(crossIdx)
    s = struct('times',[], 'count',0, 'cycleLength',[], 'rateBpm',0, 'CV_CL',NaN, ...
        'peakTimes',[], 'overshoot',[], 'MDP',[], 'amplitude',[], 'dVdtMax',[], 'APD',[], ...
        'pacemakerMechanismScore',[], 'persistentNaConductancePct',[]);
    return
end

frac  = (opts.threshold - V(crossIdx)) ./ (V(crossIdx+1) - V(crossIdx));
times = t(crossIdx) + frac .* (t(crossIdx+1) - t(crossIdx));

keepMask = true(size(times));
lastKept = 1;
for k = 2:numel(times)
    if times(k) - times(lastKept) < opts.refractory
        keepMask(k) = false;
    else
        lastKept = k;
    end
end
times    = times(keepMask);
crossIdx = crossIdx(keepMask);

nBeats = numel(times);
peakWinSamp = max(1, round(opts.peakWindow/dt));

overshoot = zeros(1,nBeats); peakTimes = zeros(1,nBeats); MDP = zeros(1,nBeats);
amplitude = zeros(1,nBeats); dVdtMax   = zeros(1,nBeats); APD = zeros(1,nBeats);
pmScore   = nan(1,nBeats);   naPct     = nan(1,nBeats);

dVdt_full = [0, diff(V)./diff(t)];

for k = 1:nBeats
    iStart = crossIdx(k);
    iEnd   = min(numel(V), iStart + peakWinSamp);
    [pk, relIdx] = max(V(iStart:iEnd));
    pkIdx = iStart + relIdx - 1;
    overshoot(k) = pk;
    peakTimes(k) = t(pkIdx);

    if k == 1
        diaStart = 1;
    else
        diaStart = find(t >= peakTimes(k-1), 1, 'first');
    end
    [mdpVal, mdpRelIdx] = min(V(diaStart:iStart));
    MDP(k)  = mdpVal;
    mdpIdx  = diaStart + mdpRelIdx - 1;
    amplitude(k) = overshoot(k) - MDP(k);

    dVdtWin = dVdt_full(max(1,iStart-2):min(numel(V),pkIdx));
    dVdtMax(k) = max(dVdtWin);

    repolV = overshoot(k) - opts.repolFraction*amplitude(k);
    searchEnd = min(numel(V), pkIdx + round(2000/dt));
    seg  = V(pkIdx:searchEnd);
    segT = t(pkIdx:searchEnd);
    fallIdx = find(seg(1:end-1) >= repolV & seg(2:end) < repolV, 1, 'first');
    if isempty(fallIdx)
        APD(k) = NaN;
    else
        f = (seg(fallIdx)-repolV) / (seg(fallIdx)-seg(fallIdx+1));
        tRepol = segT(fallIdx) + f*(segT(fallIdx+1)-segT(fallIdx));
        APD(k) = tRepol - times(k);
    end
    if k > 1 && (iStart - mdpIdx) > 2
        diaIdx = mdpIdx:iStart;
        dn = [0, diff(n(diaIdx))./diff(t(diaIdx))];
        dV = dVdt_full(diaIdx);
        consistent = (dn <= 1e-6) & (dV >= -1e-6);
        consistent = consistent(2:end); 
        pmScore(k) = 100*sum(consistent)/numel(consistent);
    end
    plWinStart = find(t >= peakTimes(k)+20, 1, 'first');
    plWinEndT  = peakTimes(k) + 480;
    if k < nBeats
        plWinEndT = min(plWinEndT, times(k+1)-20);
    end
    plWinEndT = min(plWinEndT, t(end));
    plWinEnd = find(t >= plWinEndT, 1, 'first');
    if ~isempty(plWinStart) && ~isempty(plWinEnd) && plWinEnd > plWinStart
        peakGNa = max(gNa(iStart:pkIdx));
        plateauGNa = min(gNa(plWinStart:plWinEnd));
        naPct(k) = 100*plateauGNa/peakGNa;
    end
end

s.times     = times;
s.count     = nBeats;
s.peakTimes = peakTimes;
s.overshoot = overshoot;
s.MDP       = MDP;
s.amplitude = amplitude;
s.dVdtMax   = dVdtMax;
s.APD       = APD;
s.pacemakerMechanismScore     = pmScore;
s.persistentNaConductancePct  = naPct;

if nBeats >= 2
    s.cycleLength = diff(times);
    s.rateBpm     = 60000/mean(s.cycleLength);
    s.CV_CL       = std(s.cycleLength)/mean(s.cycleLength);
else
    s.cycleLength = [];
    s.rateBpm     = 0;
    s.CV_CL       = NaN;
end
end
