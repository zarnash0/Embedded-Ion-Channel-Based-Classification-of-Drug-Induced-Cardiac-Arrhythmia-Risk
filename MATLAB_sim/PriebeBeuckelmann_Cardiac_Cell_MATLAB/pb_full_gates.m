function g = pb_full_gates(V)
xm = V+47.13;
am = zeros(size(V));
small = abs(xm) < 0.001;
am(~small) = 0.32*xm(~small) ./ (1-exp(-0.1*xm(~small)));
am(small)  = 3.2;                                  % L'Hopital limit at V=-47.13
bm = 0.08*exp(-V/11);
g.minf = am ./ (am+bm);  g.taum = 1 ./ (am+bm);

lo = V < -40;
ah = zeros(size(V)); bh = zeros(size(V));
ah(lo)  = 0.135*exp((80+V(lo))/(-6.8));
ah(~lo) = 0;
bh(lo)  = 3.56*exp(0.079*V(lo)) + 3.1e5*exp(0.35*V(lo));
bh(~lo) = 1 ./ (0.13*(1+exp(-(V(~lo)+10.66)/11.1)));
g.hinf = ah ./ (ah+bh);  g.tauh = 1 ./ (ah+bh);

aj = zeros(size(V)); bj = zeros(size(V));
aj(lo)  = ((-127140*exp(0.244*V(lo)) - 3.474e-5*exp(-0.04391*V(lo))) .* (V(lo)+37.78)) ...
          ./ (1+exp(0.311*(V(lo)+79.23)));
aj(~lo) = 0;
bj(lo)  = 0.1212*exp(-0.01052*V(lo)) ./ (1+exp(-0.1378*(V(lo)+40.14)));
bj(~lo) = 0.3*exp(-2.535e-7*V(~lo)) ./ (1+exp(-0.1*(V(~lo)+32)));
g.jinf = aj ./ (aj+bj);  g.tauj = 1 ./ (aj+bj);

ad = 14.9859*exp(-0.5*((V-22.36)/16.6813).^2) / (16.6813*sqrt(2*pi));
bd = 0.1471 - (5.3*exp(-0.5*((V-6.2744)/14.93).^2)) / (14.93*sqrt(2*pi));
g.dinf = ad ./ (ad+bd);  g.taud = 1 ./ (ad+bd);

af = 6.872e-3 ./ (1 + exp((V-6.1546)/6.1223));
bf = (0.0687*exp(-0.1081*(V+9.8255)) + 0.0112) ./ (1 + exp(-0.2779*(V+9.8255))) + 5.474e-4;
g.finf = af ./ (af+bf);  g.tauf = 1 ./ (af+bf);

ar = (0.5266*exp(-0.0166*(V-42.2912))) ./ (1 + exp(-0.0943*(V-42.2912)));
br = (5.186e-5*V + 0.5149*exp(-0.1344*(V-5.0027))) ./ (1 + exp(-0.1348*(V-5.186e-5)));
g.rinf = ar ./ (ar+br);  g.taur = 1 ./ (ar+br);

ato = (5.612e-5*V + 0.0721*exp(-0.173*(V+34.2531))) ./ (1 + exp(-0.1732*(V+34.2531)));
bto = (1.215e-4*V + 0.0767*exp(-1.66e-9*(V+34.0235))) ./ (1 + exp(-0.1604*(V+34.0235)));
g.toinf = ato ./ (ato+bto);  g.tauto = 1 ./ (ato+bto);

aXs = 3.013e-3 ./ (1 + exp((7.4454-(V+10))/14.3171));
bXs = 5.87e-3  ./ (1 + exp((5.95+V+10)/15.82));
g.Xsinf = aXs ./ (aXs+bXs);  g.tauXs = 1 ./ (aXs+bXs);

aXr = (0.005*exp(5.266e-4*(V+4.067))) ./ (1 + exp(-0.1262*(V+4.067)));
bXr = (0.016*exp(1.6e-3*(V+65.66)))  ./ (1 + exp(0.0783*(V+65.66)));
g.Xrinf = aXr ./ (aXr+bXr);  g.tauXr = 1 ./ (aXr+bXr);
end
