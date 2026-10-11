function dy = hh_derivatives(y, Ie, params)
V = y(1); m = y(2); h = y(3); n = y(4);
[a, b] = hh_gating_rates(V);

gNa_t = params.gNa * m^3 * h;
gK_t  = params.gK  * n^4;
gL_t  = params.gL;
Iion  = gL_t*(V-params.EL) + gK_t*(V-params.EK) + gNa_t*(V-params.ENa);

dV = (-Iion + Ie) / params.C;
dm = a.m*(1-m) - b.m*m;
dh = a.h*(1-h) - b.h*h;
dn = a.n*(1-n) - b.n*n;

dy = [dV; dm; dh; dn];
end
