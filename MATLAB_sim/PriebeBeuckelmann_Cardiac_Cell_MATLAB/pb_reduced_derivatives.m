function dy = pb_reduced_derivatives(y, Istim, p, const)
V=y(1); m=y(2); v=y(3); f=y(4); to=y(5); X=y(6);

g = pb_reduced_gates(V, p);
c = pb_reduced_currents(V, m, v, f, to, X, p, const, g);   

dV  = -(c.Iion - Istim);
dm  = (g.minf-m)/g.taum;
dv  = (g.vinf-v)/g.tauv;
df  = (g.finf-f)/g.tauf;
dto = (g.toinf-to)/g.tauto;
dX  = (g.Xinf-X)/g.tauX;

dy = [dV; dm; dv; df; dto; dX];
end
