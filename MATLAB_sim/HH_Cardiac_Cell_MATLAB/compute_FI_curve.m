function fi = compute_FI_curve(params, Iamps, opts)
%   Spikes are counted only after opts.settle ms following stimulus onset.
%   That alone doesn't exclude the transient: HH's near-threshold response
%   is a slowly damped oscillation that can fire past t = 120 ms and then
%   stop, so a settle window catches spikes from a cell that is no longer
%   firing. 
arguments
    params (1,1) struct
    Iamps (1,:) double
    opts.tEnd (1,1) double = 200
    opts.onset (1,1) double = 20
    opts.offset (1,1) double = 180
    opts.dt (1,1) double = 0.01
    opts.settle (1,1) double = 30
end

nI    = numel(Iamps);
rate  = zeros(1,nI);
count = zeros(1,nI);
sustained = false(1,nI);
countStart = opts.onset + opts.settle;

for k = 1:nI
    stim = make_pulse(opts.onset, opts.offset, Iamps(k));
    res  = simulate_HH(params, stim, opts.tEnd, dt=opts.dt);
    sp   = detect_spikes(res.t, res.V);
    inWindow  = sp.times >= countStart & sp.times <= opts.offset;
    count(k)  = sum(inWindow);

    w = sp.times(inWindow);
    if numel(w) >= 2
        isi = diff(w);
        rate(k) = 1000/mean(isi);
        sustained(k) = (opts.offset - w(end)) <= 2*mean(isi);
    else
        rate(k) = 0;
        sustained(k) = false;
    end
    if ~sustained(k)
        rate(k) = 0;   % not repetitive firing: report no sustained rate
    end
end

fi.I         = Iamps;
fi.rate      = rate;
fi.count     = count;
fi.sustained = sustained;
end
