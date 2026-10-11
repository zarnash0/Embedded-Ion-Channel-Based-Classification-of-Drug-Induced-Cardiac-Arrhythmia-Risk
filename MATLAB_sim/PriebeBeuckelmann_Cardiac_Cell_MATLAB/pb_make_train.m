function stimFcn = pb_make_train(period, pulseDur, amplitude, startTime)
if nargin < 4, startTime = 0; end
stimFcn = @(t) amplitude .* (mod(t-startTime, period) < pulseDur) .* (t >= startTime);
end
