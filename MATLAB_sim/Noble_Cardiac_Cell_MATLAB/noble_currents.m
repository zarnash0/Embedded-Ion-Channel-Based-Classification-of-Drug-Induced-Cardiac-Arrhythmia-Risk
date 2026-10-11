function c = noble_currents(V, m, h, n, p)
gNa = p.gNaMax .* m.^3 .* h + p.gNaBg;
gK1 = p.gK1a .* exp((-V-90)/50) + p.gK1b .* exp((V+90)/60);
gK2 = p.gK2Max .* n.^4;
gAn = p.gAn .* ones(size(V));

c.gNa = gNa;
c.gK1 = gK1;
c.gK2 = gK2;
c.gAn = gAn;

c.iNa = gNa .* (V - p.ENa);
c.iK1 = gK1 .* (V - p.EK);
c.iK2 = gK2 .* (V - p.EK);
c.iAn = gAn .* (V - p.EAn);
end
