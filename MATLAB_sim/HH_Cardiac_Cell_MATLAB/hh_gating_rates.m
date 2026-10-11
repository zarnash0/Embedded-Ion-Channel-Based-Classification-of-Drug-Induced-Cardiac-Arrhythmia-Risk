function [alpha, beta] = hh_gating_rates(V)
alpha.n = 0.01 * linoid(V, -55, 10);
beta.n  = 0.125 * exp(-(V+65)/80);

alpha.m = 0.1 * linoid(V, -40, 10);
beta.m  = 4 * exp(-(V+65)/18);

alpha.h = 0.07 * exp(-(V+65)/20);
beta.h  = 1 ./ (1 + exp(-(V+35)/10));
end

function y = linoid(V, Vhalf, k)
x = V - Vhalf;
y = zeros(size(x));
tiny = abs(x) < 1e-7;
y(~tiny) = x(~tiny) ./ (1 - exp(-x(~tiny)/k));
y(tiny)  = k;
end
