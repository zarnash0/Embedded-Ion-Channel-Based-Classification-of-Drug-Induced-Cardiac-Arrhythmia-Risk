function g = pb_reduced_gates(V, p)
xm = V + 47.13;
am = 0.32*xm ./ (1 - exp(-0.1*xm));
small = abs(xm) < 1e-6;
am(small) = 3.2;
bm = 0.08*exp(-V/11);
g.minf = am ./ (am+bm);
g.taum = 1 ./ (am+bm);

g.vinf = 0.5*(1 - tanh(7.74 + 0.12*V));
g.tauv = 0.25 + 2.24*((1 - tanh(7.74 + 0.12*V)) ./ (1 - tanh(0.07*(V+92.4))));

ad = 14.9859*exp(-0.5*((V-22.36)/16.6813).^2) / (16.6813*sqrt(2*pi));
bd = 0.1471 - (5.3*exp(-0.5*((V-6.2744)/14.93).^2)) / (14.93*sqrt(2*pi));
g.dinf = ad ./ (ad+bd);
af = 6.872e-3 ./ (1 + exp((V-6.1546)/6.1223));
bf = (0.0687*exp(-0.1081*(V+9.8255)) + 0.0112) ./ (1 + exp(-0.2779*(V+9.8255))) + 5.474e-4;
g.finf = af ./ (af+bf);
g.tauf = 1 ./ (af+bf);

ar = (0.5266*exp(-0.0166*(V-42.2912))) ./ (1 + exp(-0.0943*(V-42.2912)));
br = (5.186e-5*V + 0.5149*exp(-0.1344*(V-5.0027))) ./ (1 + exp(-0.1348*(V-5.186e-5)));
g.rinf = ar ./ (ar+br);

Vs = V - p.Vshift;
ato_base = (5.612e-5*Vs + 0.0721*exp(-0.173*(Vs+34.2531))) ./ (1 + exp(-0.1732*(Vs+34.2531)));
bto_base = (1.215e-4*Vs + 0.0767*exp(-1.66e-9*(Vs+34.0235))) ./ (1 + exp(-0.1604*(Vs+34.0235)));
g.toinf = ato_base ./ (ato_base+bto_base);          
g.tauto = p.invP ./ (ato_base+bto_base);            

if p.cellType == "M"
    g.Xinf = 0.972 ./ (1 + exp(-2.036 - 0.0834*V));
    g.tauX = 380*exp(-((25.5+V).^2)/156) + 166*(1+tanh(0.558+0.0169*V));
else
    g.Xinf = 0.988 ./ (1 + exp(-0.861 - 0.0620*V));

    tauXprime = 40*(1 - tanh(160 + 2*V));
    g.tauX = 240*exp(-((25.5+V).^2)/156) + 182*(1+tanh(0.154+0.0116*V)) + tauXprime;
end

end
