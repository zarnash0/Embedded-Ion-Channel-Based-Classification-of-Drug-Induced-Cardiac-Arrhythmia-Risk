function stimFcn = make_pulse(onset, offset, amplitude)
%   pulse of amplitude (nA/mm^2) active for ONSET <= t < OFFSET (ms).
stimFcn = @(t) amplitude .* (t >= onset & t < offset);
end
