function dy = pb_full_derivatives(y, Istim, p)

V=y(1); m=y(2); h=y(3); j=y(4); d=y(5); f=y(6); r=y(7); to=y(8);
Xs=y(9); Xr=y(10); Nai=y(11); Ki=y(12); Cai=y(13); CaJSR=y(14); CaNSR=y(15);
APtrack=y(16); APtrack2=y(17); Cainfluxtrack=y(18);

g = pb_full_gates(V);
c = pb_full_currents(y(1:15), p);

dV  = -c.Iion + Istim;
dm  = (g.minf-m)/g.taum;
dh  = (g.hinf-h)/g.tauh;
dj  = (g.jinf-j)/g.tauj;
dd  = (g.dinf-d)/g.taud;
df  = (g.finf-f)/g.tauf;
dr  = (g.rinf-r)/g.taur;
dto = (g.toinf-to)/g.tauto;
dXs = (g.Xsinf-Xs)/g.tauXs;
dXr = (g.Xrinf-Xr)/g.tauXr;

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
    dCainfluxtrack = -p.Cm*((c.iCa-c.iNaCa)+c.iCab) / (2*p.Vmyo*p.F);
elseif (APtrack2 > 0.01) && (APtrack <= 0.2)
    dCainfluxtrack = 0;
else
    dCainfluxtrack = -0.5*Cainfluxtrack;
end

if Cainfluxtrack > p.delta_Ca_ith
    Grel = p.Grel_max*(Cainfluxtrack-p.delta_Ca_ith)/(p.Kmrel+Cainfluxtrack-p.delta_Ca_ith) * (1-APtrack2)*APtrack2;
else
    Grel = 0;
end
irel = Grel*(CaJSR-Cai);

Cai_bufc   = 1/(1+p.buffon*((p.CMDNmax*p.KmCMDN/(p.KmCMDN+Cai)^2) + (p.Tnmax*p.KmTn/(p.KmTn+Cai)^2)));
CaJSR_bufc = 1/(1+(p.CSQNmax*p.KmCSQN/(p.KmCSQN+CaJSR)^2));

% NOTE ON SIGN: the CellML's Ki balance sums i_Stim using the SAME
% outward-positive (LR) convention as every other current; this file's
% Istim is defined the opposite way (positive = depolarizing), so the
% substitution here is i_Stim_LR = -Istim.
dNai   = -p.Cm*(c.iNa + c.iNab + 3*c.iNaCa + 3*c.iNaK) / (p.Vmyo*p.F);
dKi    = -p.Cm*((c.ito + c.iKr + c.iK1 + (-Istim) + c.iKs) - 2*c.iNaK) / (p.Vmyo*p.F);
dCai   = Cai_bufc*( -p.Cm*((c.iCa-2*c.iNaCa)+c.iCab)/(2*p.Vmyo*p.F) + irel*p.VJSR/p.Vmyo + (c.ileak-c.iup)*p.VNSR/p.Vmyo );
dCaJSR = CaJSR_bufc*(c.itr - irel);
dCaNSR = -(c.ileak + (p.VJSR/p.VNSR)*c.itr - c.iup);

dy = [dV;dm;dh;dj;dd;df;dr;dto;dXs;dXr;dNai;dKi;dCai;dCaJSR;dCaNSR;dAPtrack;dAPtrack2;dCainfluxtrack];
end
