function stimFcn = noble_make_pulse(onset, offset, amplitude)
stimFcn = @(t) amplitude .* (t >= onset & t < offset);
end
