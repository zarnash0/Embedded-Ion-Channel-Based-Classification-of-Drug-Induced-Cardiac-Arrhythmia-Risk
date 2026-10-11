function dy = noble_derivatives(y, Istim, p)
V = y(1); m = y(2); h = y(3); n = y(4);
[a, b] = noble_gating_rates(V);
c = noble_currents(V, m, h, n, p);

dV = -(c.iNa + c.iK1 + c.iK2 + c.iAn - Istim) / p.Cm;   % mV/s
dm = a.m*(1-m) - b.m*m;                                  % 1/s
dh = a.h*(1-h) - b.h*h;                                  % 1/s
dn = a.n*(1-n) - b.n*n;                                  % 1/s

dy = [dV; dm; dh; dn];
end
