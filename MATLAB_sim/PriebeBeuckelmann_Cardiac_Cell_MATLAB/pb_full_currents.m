function c = pb_full_currents(y, p)
V=y(1); m=y(2); h=y(3); j=y(4); d=y(5); f=y(6); r=y(7); to=y(8);
Xs=y(9); Xr=y(10); Nai=y(11); Ki=y(12); Cai=y(13); CaJSR=y(14); CaNSR=y(15);

RTF = p.RTF;
ENa = RTF*log(p.Nao/Nai);
ECa = RTF/2*log(p.Cao/Cai);
Eto = RTF*log((0.043*p.Nao+p.Ko)/(0.043*Nai+Ki));
EKs = RTF*log((0.01833*p.Nao+p.Ko)/(0.01833*Nai+Ki));
EK  = RTF*log(p.Ko/Ki);

fCa = p.KmCa_fCa/(p.KmCa_fCa+Cai);

c.iNa = p.gNa*m^3*h*j*(V-ENa);
c.iCa = p.gCa*d*f*fCa*(V-ECa);
c.ito = p.gto*r*to*(V-Eto);
c.iKs = p.gKs*Xs^2*(V-EKs);

rik = 1/(1+exp((V+26)/23));
c.iKr = p.gKr*Xr*rik*(V-EK);

aK1 = 0.1/(1+exp(0.06*(V-EK-200)));
bK1 = (3*exp(2e-4*(V-EK+100)) + exp(0.1*(V-EK-10))) / (1+exp(-0.5*(V-EK)));
K1inf = aK1/(aK1+bK1);
c.iK1 = p.gK1*K1inf*(V-EK);

c.iNab = p.gNab*(V-ENa);
c.iCab = p.gCab*(V-ECa);

sigma = (1/7)*(exp(p.Nao/67.3)-1);
fNaK = 1/(1+0.1245*exp(-0.1*V/RTF)+0.0365*sigma*exp(-V/RTF));
c.iNaK = p.INaKmax*fNaK*(1/(1+(p.KmNai_NaK/Nai)^1.5))*(p.Ko/(p.Ko+p.KmKo_NaK));

fNaCa = (1/((p.KmNa_NaCa^3+p.Nao^3)*(p.KmCa_NaCa+p.Cao)*(1+p.ksat_NaCa*exp((p.eta_NaCa-1)*V/RTF)))) ...
        * (exp(p.eta_NaCa*V/RTF)*Nai^3*p.Cao - exp((p.eta_NaCa-1)*V/RTF)*p.Nao^3*Cai);
c.iNaCa = p.KNaCa*fNaCa;

c.iup   = p.Iupmax*Cai/(Cai+p.Kmup);
c.ileak = p.Kleak*CaNSR;
c.itr   = (CaNSR-CaJSR)/p.tautr;

c.ENa=ENa; c.ECa=ECa; c.Eto=Eto; c.EKs=EKs; c.EK=EK;

c.Iion = c.iNa + c.iCa + c.ito + c.iKr + c.iKs + c.iK1 + c.iNaK + c.iNaCa + c.iNab + c.iCab;
end
