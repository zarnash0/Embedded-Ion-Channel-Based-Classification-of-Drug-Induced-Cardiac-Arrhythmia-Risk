function out = simulate_Noble(params, tEnd, opts)
%   INTEGRATION SCHEME: at each step, freezes m, h, n (as usual) and
%   freezes the purely voltage-dependent i_K1 conductance at the current
%   V(i)- exactly the same "frozen conductance, exact-solve the
%   resulting locally-affine ODE" treatment already applied to the
%   gated conductances -- so the whole current sum becomes affine in V
%   for the step, and V, m, h, n all update via
%   x(t+dt) = xinf + (x(t)-xinf)*exp(-dt/tau(x)), unconditionally stable
%   for any dt>0 (though accuracy of the fast upstroke still needs a
%   sufficiently small dt -- see the dt/10 convergence check in
%   run_Noble_demo.m). "rk4" integrates the raw ODE system directly
%   as an independent correctness cross-check.

arguments
    params (1,1) struct
    tEnd (1,1) double {mustBePositive}
    opts.dt (1,1) double {mustBePositive} = 0.005
    opts.stim = []
    opts.V0 (1,1) double = params.V0
    opts.m0 (1,1) double = params.m0
    opts.h0 (1,1) double = params.h0
    opts.n0 (1,1) double = params.n0
    opts.method (1,1) string {mustBeMember(opts.method, ["exponential","rk4"])} = "exponential"
    opts.returnCurrents (1,1) logical = true
end

t  = 0:opts.dt:tEnd;
N  = numel(t);
dt_s = opts.dt/1000;

if isempty(opts.stim)
    Istim = zeros(1,N);
elseif isa(opts.stim, 'function_handle')
    Istim = opts.stim(t);
else
    Istim = opts.stim(:)';
    if numel(Istim) ~= N
        error('simulate_Noble:stimLength', ...
            'Numeric stim must have %d samples to match 0:dt:tEnd (got %d).', N, numel(Istim));
    end
end

V = zeros(1,N); m = zeros(1,N); h = zeros(1,N); n = zeros(1,N);
if opts.returnCurrents
    iNa = zeros(1,N); iK1 = zeros(1,N); iK2 = zeros(1,N); iAn = zeros(1,N);
    gNa = zeros(1,N); gK1 = zeros(1,N); gK2 = zeros(1,N); gAn = zeros(1,N);
else
    iNa = nan(1,N);   iK1 = nan(1,N);   iK2 = nan(1,N);   iAn = nan(1,N);
    gNa = nan(1,N);   gK1 = nan(1,N);   gK2 = nan(1,N);   gAn = nan(1,N);
end

V(1) = opts.V0; m(1) = opts.m0; h(1) = opts.h0; n(1) = opts.n0;
keepI = opts.returnCurrents;

if opts.method == "exponential"
    for i = 1:N-1
        [a, b] = noble_gating_rates(V(i));
        minf = a.m/(a.m+b.m); taum_s = 1/(a.m+b.m);
        hinf = a.h/(a.h+b.h); tauh_s = 1/(a.h+b.h);
        ninf = a.n/(a.n+b.n); taun_s = 1/(a.n+b.n);

        cur = noble_currents(V(i), m(i), h(i), n(i), params);
        if keepI
            iNa(i) = cur.iNa; iK1(i) = cur.iK1; iK2(i) = cur.iK2; iAn(i) = cur.iAn;
            gNa(i) = cur.gNa; gK1(i) = cur.gK1; gK2(i) = cur.gK2; gAn(i) = cur.gAn;
        end

        gTot = cur.gNa + cur.gK1 + cur.gK2 + cur.gAn;
        Vinf = (cur.gNa*params.ENa + cur.gK1*params.EK + cur.gK2*params.EK + cur.gAn*params.EAn + Istim(i)) / gTot;
        tauV_s = params.Cm / gTot;

        m(i+1) = minf + (m(i)-minf)*exp(-dt_s/taum_s);
        h(i+1) = hinf + (h(i)-hinf)*exp(-dt_s/tauh_s);
        n(i+1) = ninf + (n(i)-ninf)*exp(-dt_s/taun_s);
        V(i+1) = Vinf + (V(i)-Vinf)*exp(-dt_s/tauV_s);
    end
    if keepI
        cur = noble_currents(V(N), m(N), h(N), n(N), params);
        iNa(N) = cur.iNa; iK1(N) = cur.iK1; iK2(N) = cur.iK2; iAn(N) = cur.iAn;
        gNa(N) = cur.gNa; gK1(N) = cur.gK1; gK2(N) = cur.gK2; gAn(N) = cur.gAn;
    end
else
    for i = 1:N-1
        y = [V(i); m(i); h(i); n(i)];
        k1 = noble_derivatives(y,             Istim(i), params);
        k2 = noble_derivatives(y+0.5*dt_s*k1, Istim(i), params);
        k3 = noble_derivatives(y+0.5*dt_s*k2, Istim(i), params);
        k4 = noble_derivatives(y+dt_s*k3,     Istim(i), params);
        ynext = y + (dt_s/6)*(k1 + 2*k2 + 2*k3 + k4);
        V(i+1) = ynext(1); m(i+1) = ynext(2); h(i+1) = ynext(3); n(i+1) = ynext(4);
    end
    if keepI
        for i = 1:N
            cur = noble_currents(V(i), m(i), h(i), n(i), params);
            iNa(i) = cur.iNa; iK1(i) = cur.iK1; iK2(i) = cur.iK2; iAn(i) = cur.iAn;
            gNa(i) = cur.gNa; gK1(i) = cur.gK1; gK2(i) = cur.gK2; gAn(i) = cur.gAn;
        end
    end
end

out = struct('t',t, 'V',V, 'm',m, 'h',h, 'n',n, 'Istim',Istim, ...
             'iNa',iNa, 'iK1',iK1, 'iK2',iK2, 'iAn',iAn, ...
             'gNa',gNa, 'gK1',gK1, 'gK2',gK2, 'gAn',gAn, ...
             'params',params, 'dt',opts.dt, 'method',char(opts.method));
end
