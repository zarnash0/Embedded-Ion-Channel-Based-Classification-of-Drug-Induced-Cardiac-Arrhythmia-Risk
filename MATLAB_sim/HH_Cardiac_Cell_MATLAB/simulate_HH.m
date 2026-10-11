function out = simulate_HH(params, stim, tEnd, opts)
%   The default "exponential" method applies, at every step and to V, m,
%   h, and n alike, the exact solution of the locally-linearized ODE
%   x(t+dt) = xinf + (x(t)-xinf)*exp(-dt/tau(x)) -- the same scheme used
%   in the Cornell BIONB 2220 Hodgkin-Huxley lab (B. Land) that the
%   hardware cardiac-cell project cites as its MATLAB starting point
%   This update is unconditionally stable for any dt>0, unlike forward
%   Euler. "rk4" is provided purely as an independent cross-check of the
%   equations/parameters using classical 4th-order Runge-Kutta on the raw
%   ODE system.

arguments
    params (1,1) struct
    stim
    tEnd (1,1) double {mustBePositive}
    opts.dt (1,1) double {mustBePositive} = 0.01
    opts.V0 (1,1) double = params.Vrest
    opts.method (1,1) string {mustBeMember(opts.method, ["exponential","rk4"])} = "exponential"
    opts.returnCurrents (1,1) logical = true
end

t = 0:opts.dt:tEnd;
N = numel(t);

if isa(stim, 'function_handle')
    Ie = stim(t);
else
    Ie = stim(:)';
    if numel(Ie) ~= N
        error('simulate_HH:stimLength', ...
            'Numeric stim must have %d samples to match 0:dt:tEnd (got %d).', N, numel(Ie));
    end
end

V = zeros(1,N); m = zeros(1,N); h = zeros(1,N); n = zeros(1,N);
if opts.returnCurrents
    INa = zeros(1,N); IK = zeros(1,N); IL = zeros(1,N);
else
    INa = nan(1,N);   IK = nan(1,N);   IL = nan(1,N);
end

V(1) = opts.V0;
[a0, b0] = hh_gating_rates(V(1));
m(1) = a0.m/(a0.m+b0.m);
h(1) = a0.h/(a0.h+b0.h);
n(1) = a0.n/(a0.n+b0.n);

dt = opts.dt;
keepI = opts.returnCurrents;

if opts.method == "exponential"
    for i = 1:N-1
        [a, b] = hh_gating_rates(V(i));
        minf = a.m/(a.m+b.m); taum = 1/(a.m+b.m);
        hinf = a.h/(a.h+b.h); tauh = 1/(a.h+b.h);
        ninf = a.n/(a.n+b.n); taun = 1/(a.n+b.n);

        gNa_t = params.gNa * m(i)^3 * h(i);
        gK_t  = params.gK  * n(i)^4;
        gL_t  = params.gL;

        if keepI
            INa(i) = gNa_t*(V(i)-params.ENa);
            IK(i)  = gK_t *(V(i)-params.EK);
            IL(i)  = gL_t *(V(i)-params.EL);
        end

        gTot = gL_t + gK_t + gNa_t;
        Vinf = (gL_t*params.EL + gK_t*params.EK + gNa_t*params.ENa + Ie(i)) / gTot;
        tauV = params.C / gTot;

        m(i+1) = minf + (m(i)-minf)*exp(-dt/taum);
        h(i+1) = hinf + (h(i)-hinf)*exp(-dt/tauh);
        n(i+1) = ninf + (n(i)-ninf)*exp(-dt/taun);
        V(i+1) = Vinf + (V(i)-Vinf)*exp(-dt/tauV);
    end
    if keepI
        gNa_t = params.gNa*m(N)^3*h(N); gK_t = params.gK*n(N)^4;
        INa(N) = gNa_t*(V(N)-params.ENa);
        IK(N)  = gK_t *(V(N)-params.EK);
        IL(N)  = params.gL*(V(N)-params.EL);
    end
else
    for i = 1:N-1
        y = [V(i); m(i); h(i); n(i)];
        k1 = hh_derivatives(y,              Ie(i), params);
        k2 = hh_derivatives(y + 0.5*dt*k1,  Ie(i), params);
        k3 = hh_derivatives(y + 0.5*dt*k2,  Ie(i), params);
        k4 = hh_derivatives(y + dt*k3,      Ie(i), params);
        ynext = y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        V(i+1) = ynext(1); m(i+1) = ynext(2); h(i+1) = ynext(3); n(i+1) = ynext(4);
    end
    if keepI
        for i = 1:N
            gNa_t = params.gNa*m(i)^3*h(i); gK_t = params.gK*n(i)^4;
            INa(i) = gNa_t*(V(i)-params.ENa);
            IK(i)  = gK_t *(V(i)-params.EK);
            IL(i)  = params.gL*(V(i)-params.EL);
        end
    end
end

out = struct('t',t, 'V',V, 'm',m, 'h',h, 'n',n, 'Ie',Ie, ...
             'INa',INa, 'IK',IK, 'IL',IL, ...
             'params',params, 'dt',dt, 'method',char(opts.method));
end
