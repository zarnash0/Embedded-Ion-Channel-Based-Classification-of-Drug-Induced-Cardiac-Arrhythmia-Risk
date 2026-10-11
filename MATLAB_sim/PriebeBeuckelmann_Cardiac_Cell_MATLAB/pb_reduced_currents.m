function c = pb_reduced_currents(V, m, v, f, to, X, p, const, g)
RTF = const.RTF;

ENa = RTF*log(p.Nae/p.Nai);
ECa = RTF/2*log(p.Cae/p.Cai);
Eto = RTF*log((0.043*p.Nae+p.Ke)/(0.043*p.Nai+p.Ki));
EK  = RTF*log(p.Ke/p.Ki);

if nargin < 9 || isempty(g)
    g = pb_reduced_gates(V, p);
end

c.INa = p.gNa * m.^3 .* v.^2 .* (V-ENa);
c.ICa = p.gCa * g.dinf .* f .* p.fCa .* (V-ECa);
c.Ito = p.gto * g.rinf .* to .* (V-Eto);
c.IK  = p.gK  * X.^2 .* (V-EK);

aK1 = 0.1 ./ (1 + exp(0.06*(V-EK-200)));
bK1 = (3*exp(2e-4*(V-EK+100)) + exp(0.1*(V-EK-10))) ./ (1 + exp(-0.5*(V-EK)));
K1inf = aK1 ./ (aK1+bK1);
c.IK1 = p.gK1 * K1inf .* (V-EK);

c.INab = p.gNab * (V-ENa);
c.ICab = p.gCab * (V-ECa);

sigma = (1/7)*(exp(p.Nae/67.3)-1);
fNaK  = 1 ./ (1 + 0.1245*exp(-0.1*V/RTF) + 0.0365*sigma.*exp(-V/RTF));
fNaK2 = 1 ./ (1 + (p.KmNai_NaK/p.Nai)^1.5) * (p.Ke/(p.Ke+p.KmKo_NaK));
c.INaK = p.gNaK * fNaK * fNaK2;

fNaCa = (p.KmNa_NaCa^3+p.Nae^3)^-1 * (p.KmCa_NaCa+p.Cae)^-1 ...
    ./ (1 + p.ksat_NaCa*exp((p.eta_NaCa-1)*V/RTF)) ...
    .* ( exp(p.eta_NaCa*V/RTF)*p.Nai^3*p.Cae - exp((p.eta_NaCa-1)*V/RTF)*p.Nae^3*p.Cai );
c.INaCa = p.gNaCa * fNaCa;

c.Iion = c.INa + c.ICa + c.Ito + c.IK + c.IK1 + c.INab + c.ICab + c.INaK + c.INaCa;

c.ENa = ENa; c.ECa = ECa; c.Eto = Eto; c.EK = EK;
end
