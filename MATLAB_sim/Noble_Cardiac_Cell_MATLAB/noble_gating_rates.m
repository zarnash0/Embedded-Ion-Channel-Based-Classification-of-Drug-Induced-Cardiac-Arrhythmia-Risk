function [alpha, beta] = noble_gating_rates(V)
alpha.m = 100 * ratExp(-V-48, 15);
beta.m  = 120 * ratExp(V+8, 5);

alpha.h = 170 * exp((-V-90)/20);
beta.h  = 1000 ./ (1 + exp((-V-42)/10));

alpha.n = 0.1 * ratExp(-V-50, 10);
beta.n  = 2 * exp((-V-90)/80);
end

function y = ratExp(x, k)
y = zeros(size(x));
tiny = abs(x) < 1e-7;
y(~tiny) = x(~tiny) ./ (exp(x(~tiny)/k) - 1);
y(tiny)  = k;
end
