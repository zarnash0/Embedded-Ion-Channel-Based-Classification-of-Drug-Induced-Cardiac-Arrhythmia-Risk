function s = detect_spikes(t, V, opts)
arguments
    t (1,:) double
    V (1,:) double
    opts.threshold (1,1) double = -10
    opts.refractory (1,1) double = 1
    opts.window (1,1) double = 15
end

crossIdx = find(V(1:end-1) < opts.threshold & V(2:end) >= opts.threshold);

if isempty(crossIdx)
    s = struct('times',[], 'count',0, 'ISI',[], 'rate',0, 'CV_ISI',NaN, ...
                'peakTimes',[], 'amplitude',[], 'halfWidth',[]);
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

dt        = t(2) - t(1);
winSamp   = max(1, round(opts.window/dt));
baseline  = V(1);
nSpk      = numel(times);
amplitude = zeros(1,nSpk);
halfWidth = zeros(1,nSpk);
peakTimes = zeros(1,nSpk);

for k = 1:nSpk
    iStart = crossIdx(k);
    iEnd = min(numel(V), iStart + winSamp);
    if k < nSpk
        iEnd = min(iEnd, crossIdx(k+1) - 1);
    end
    iEnd = max(iEnd, iStart);
    [pk, relIdx] = max(V(iStart:iEnd));
    pkIdx = iStart + relIdx - 1;
    amplitude(k) = pk - baseline;
    peakTimes(k) = t(pkIdx);

    halfV = baseline + amplitude(k)/2;

    riseSeg = V(max(1,iStart-winSamp):pkIdx);
    riseT   = t(max(1,iStart-winSamp):pkIdx);
    riseIdx = find(riseSeg(1:end-1) < halfV & riseSeg(2:end) >= halfV, 1, 'last');
    if isempty(riseIdx)
        tRise = NaN;   
    else
        f = (halfV-riseSeg(riseIdx))/(riseSeg(riseIdx+1)-riseSeg(riseIdx));
        tRise = riseT(riseIdx) + f*(riseT(riseIdx+1)-riseT(riseIdx));
    end

    fallEnd  = min(numel(V), pkIdx+winSamp);
    fallSeg  = V(pkIdx:fallEnd);
    fallT    = t(pkIdx:fallEnd);
    fallIdx  = find(fallSeg(1:end-1) >= halfV & fallSeg(2:end) < halfV, 1, 'first');
    if isempty(fallIdx)
        tFall = NaN;
    else
        f = (fallSeg(fallIdx)-halfV)/(fallSeg(fallIdx)-fallSeg(fallIdx+1));
        tFall = fallT(fallIdx) + f*(fallT(fallIdx+1)-fallT(fallIdx));
    end

    halfWidth(k) = tFall - tRise;
end

s.times     = times;
s.count     = nSpk;
s.peakTimes = peakTimes;
s.amplitude = amplitude;
s.halfWidth = halfWidth;

if nSpk >= 2
    s.ISI    = diff(times);
    s.rate   = 1000/mean(s.ISI);
    s.CV_ISI = std(s.ISI)/mean(s.ISI);
else
    s.ISI    = [];
    s.rate   = 0;
    s.CV_ISI = NaN;
end
end
