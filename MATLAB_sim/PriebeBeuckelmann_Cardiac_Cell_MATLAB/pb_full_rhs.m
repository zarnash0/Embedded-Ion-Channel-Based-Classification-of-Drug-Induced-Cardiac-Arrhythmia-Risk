function dy = pb_full_rhs(y, Istim, p)
V=y(1); m=y(2); h=y(3); jg=y(4); d=y(5); f=y(6); r=y(7); to=y(8);
Xs=y(9); Xr=y(10); Nai=y(11); Ki=y(12); Cai=y(13); CaJSR=y(14); CaNSR=y(15);
APtrack=y(16); APtrack2=y(17); Cainfluxtrack=y(18);

RTF = p.RTF;

%  reversal potentials 
ENa = RTF*log(p.Nao/Nai);
ECa = 0.5*RTF*log(p.Cao/Cai);
Eto = RTF*log(p.EtoNum/(0.043*Nai+Ki));
EKs = RTF*log(p.EKsNum/(0.01833*Nai+Ki));
EK  = RTF*log(p.Ko/Ki);

%  gates 
xm = V+47.13;
if abs(xm) < 0.001
    am = 3.2;
else
    am = 0.32*xm/(1-exp(-0.1*xm));
end
bm = 0.08*exp(-V/11);

if V < -40
    ah = 0.135*exp((80+V)/(-6.8));
    bh = 3.56*exp(0.079*V) + 3.1e5*exp(0.35*V);
    aj = ((-127140*exp(0.244*V) - 3.474e-5*exp(-0.04391*V))*(V+37.78)) / (1+exp(0.311*(V+79.23)));
    bj = 0.1212*exp(-0.01052*V) / (1+exp(-0.1378*(V+40.14)));
else
    ah = 0;
    bh = 1/(0.13*(1+exp(-(V+10.66)/11.1)));
    aj = 0;
    bj = 0.3*exp(-2.535e-7*V) / (1+exp(-0.1*(V+32)));
end

ad = p.dCoefA*exp(-0.5*((V-22.36)/16.6813)^2);
bd = 0.1471 - p.dCoefB*exp(-0.5*((V-6.2744)/14.93)^2);

af = 6.872e-3 / (1 + exp((V-6.1546)/6.1223));
bf = (0.0687*exp(-0.1081*(V+9.8255)) + 0.0112) / (1 + exp(-0.2779*(V+9.8255))) + 5.474e-4;

ar = (0.5266*exp(-0.0166*(V-42.2912))) / (1 + exp(-0.0943*(V-42.2912)));
br = (5.186e-5*V + 0.5149*exp(-0.1344*(V-5.0027))) / (1 + exp(-0.1348*(V-5.186e-5)));

ato = (5.612e-5*V + 0.0721*exp(-0.173*(V+34.2531))) / (1 + exp(-0.1732*(V+34.2531)));
bto = (1.215e-4*V + 0.0767*exp(-1.66e-9*(V+34.0235))) / (1 + exp(-0.1604*(V+34.0235)));

aXs = 3.013e-3 / (1 + exp((7.4454-(V+10))/14.3171));
bXs = 5.87e-3  / (1 + exp((5.95+V+10)/15.82));

aXr = (0.005*exp(5.266e-4*(V+4.067))) / (1 + exp(-0.1262*(V+4.067)));
bXr = (0.016*exp(1.6e-3*(V+65.66)))  / (1 + exp(0.0783*(V+65.66)));

% currents 
fCa = p.KmCa_fCa/(p.KmCa_fCa+Cai);

iNa = p.gNa*m*m*m*h*jg*(V-ENa);
iCa = p.gCa*d*f*fCa*(V-ECa);
ito = p.gto*r*to*(V-Eto);
iKs = p.gKs*Xs*Xs*(V-EKs);

rik = 1/(1+exp((V+26)/23));
iKr = p.gKr*Xr*rik*(V-EK);

aK1 = 0.1/(1+exp(0.06*(V-EK-200)));
bK1 = (3*exp(2e-4*(V-EK+100)) + exp(0.1*(V-EK-10))) / (1+exp(-0.5*(V-EK)));
iK1 = p.gK1*(aK1/(aK1+bK1))*(V-EK);

iNab = p.gNab*(V-ENa);
iCab = p.gCab*(V-ECa);

fNaK = 1/(1+0.1245*exp(-0.1*V/RTF)+p.sigmaNaK0365*exp(-V/RTF));
iNaK = p.INaKmax*fNaK*(1/(1+(p.KmNai_NaK/Nai)^1.5))*p.KoFrac;

iNaCa = p.KNaCa*(1/(p.NaCaDen*(1+p.ksat_NaCa*exp((p.eta_NaCa-1)*V/RTF)))) ...
        * (exp(p.eta_NaCa*V/RTF)*Nai*Nai*Nai*p.Cao - exp((p.eta_NaCa-1)*V/RTF)*p.Nao3*Cai);

iup   = p.Iupmax*Cai/(Cai+p.Kmup);
ileak = p.Kleak*CaNSR;
itr   = (CaNSR-CaJSR)/p.tautr;

Iion = iNa + iCa + ito + iKr + iKs + iK1 + iNaK + iNaCa + iNab + iCab;

% CICR tracking 
dV = -Iion + Istim;

if dV > 150
    dAPtrack = 100*(1-APtrack) - 0.5*APtrack;
else
    dAPtrack = -0.5*APtrack;
end

if (APtrack < 0.2) && (APtrack > 0.18)
    dAPtrack2 = 100*(1-APtrack2) - 0.5*APtrack2;
else
    dAPtrack2 = -0.5*APtrack2;
end

if APtrack > 0.2
    dCainfluxtrack = -p.Cm*((iCa-iNaCa)+iCab) / (2*p.Vmyo*p.F);
elseif (APtrack2 > 0.01) && (APtrack <= 0.2)
    dCainfluxtrack = 0;
else
    dCainfluxtrack = -0.5*Cainfluxtrack;
end

if Cainfluxtrack > p.delta_Ca_ith
    Grel = p.Grel_max*(Cainfluxtrack-p.delta_Ca_ith)/(p.Kmrel+Cainfluxtrack-p.delta_Ca_ith)*(1-APtrack2)*APtrack2;
else
    Grel = 0;
end
irel = Grel*(CaJSR-Cai);

Cai_bufc   = 1/(1+p.buffon*((p.CMDNmax*p.KmCMDN/(p.KmCMDN+Cai)^2) + (p.Tnmax*p.KmTn/(p.KmTn+Cai)^2)));
CaJSR_bufc = 1/(1+(p.CSQNmax*p.KmCSQN/(p.KmCSQN+CaJSR)^2));

%
dy = zeros(18,1);
dy(1)  = dV;
dy(2)  = am*(1-m) - bm*m;
dy(3)  = ah*(1-h) - bh*h;
dy(4)  = aj*(1-jg) - bj*jg;
dy(5)  = ad*(1-d) - bd*d;
dy(6)  = af*(1-f) - bf*f;
dy(7)  = ar*(1-r) - br*r;
dy(8)  = ato*(1-to) - bto*to;
dy(9)  = aXs*(1-Xs) - bXs*Xs;
dy(10) = aXr*(1-Xr) - bXr*Xr;
dy(11) = -p.Cm*(iNa + iNab + 3*iNaCa + 3*iNaK) / (p.Vmyo*p.F);
dy(12) = -p.Cm*((ito + iKr + iK1 + (-Istim) + iKs) - 2*iNaK) / (p.Vmyo*p.F);
dy(13) = Cai_bufc*( -p.Cm*((iCa-2*iNaCa)+iCab)/(2*p.Vmyo*p.F) + irel*p.VJSR/p.Vmyo + (ileak-iup)*p.VNSR/p.Vmyo );
dy(14) = CaJSR_bufc*(itr - irel);
dy(15) = -(ileak + (p.VJSR/p.VNSR)*itr - iup);
dy(16) = dAPtrack;
dy(17) = dAPtrack2;
dy(18) = dCainfluxtrack;
end
