function out = simulate_PB_reduced(params, tEnd, opts)
%SIMULATE_PB_REDUCED Integrate the Bernus et al. (2002) 6-variable reduced model.
%   OUT = SIMULATE_PB_REDUCED(PARAMS, TEND) integrates from t=0 to
%   t=TEND (ms), paced by a periodic stimulus train (default: 1 Hz,
%   3-ms pulses of 15 pA/pF, matching the Priebe-Beuckelmann CellML
%   reference simulation's stimulus protocol -- note the CellML's own
%   sign convention has i_Stim NEGATIVE for a depolarizing/inward
%   stimulus, since it is summed directly into I_ion with the same
%   outward-positive convention as every other current; this code instead
%   treats Istim as positive-depolarizing throughout, entering via +Istim
%   in the membrane equation, so the magnitude here (+15) is the same
%   stimulus with the sign flipped to match this file's own convention).
%
%   Name-value options:
%     dt      time step in ms                                [default 0.02]
%     stim    function handle stim(t)->pA/pF, or a numeric     [default 1 Hz train]
%             vector already sampled on 0:dt:tEnd
%     V0      initial membrane potential                 [default params.V0]
%     method  "exponential" (default) or "rk4"
%
%   INTEGRATION SCHEME ("exponential", default): the five gates (m, v, f,
%   to, X) are all classical Hodgkin-Huxley-type relaxation equations and
%   are updated via Rush-Larsen's exponential update x(t+dt) =
%   xinf+(x-xinf)*exp(-dt/tau) -- explicitly the technique both source
%   papers cite for gating-variable integration (Rush S & Larsen H,
%   1978). The membrane potential is updated with the same scheme after
%   freezing, at each step, (a) the gated/instantaneous conductances at
%   their current values (as usual) and (b) the two currents that are
%   not of simple conductance*(V-E) form -- INaK and INaCa, both genuine
%   nonlinear functions of V -- as constants for that step (the same
%   treatment as an external stimulus). This is a standard operator-split
%   semi-implicit scheme; "rk4" integrates the raw ODE system directly
%   (via PB_REDUCED_DERIVATIVES) as an independent correctness check.

arguments
    params (1,1) struct
    tEnd (1,1) double {mustBePositive}
    opts.dt (1,1) double {mustBePositive} = 0.02
    opts.stim = []
    opts.V0 (1,1) double = params.V0
    opts.method (1,1) string {mustBeMember(opts.method, ["exponential","rk4"])} = "exponential"
    opts.returnCurrents (1,1) logical = true
end
% returnCurrents exists so COMPARE_PB_MODELS can time both models with
% output bookkeeping switched off; see SIMULATE_PB_FULL. In this model
% the currents fall out of the integration loop almost for free, so the
% flag changes little here -- which is exactly why it must also be
% honoured here, otherwise the comparison would remove the cost from one
% model only.

const = pb_constants();
t = 0:opts.dt:tEnd;
N = numel(t);
dt = opts.dt;

if isempty(opts.stim)
    stimFcn = pb_make_train(1000, 3, 15, 50);   % 50 ms offset, see SIMULATE_PB_FULL
    Istim = stimFcn(t);
elseif isa(opts.stim, 'function_handle')
    Istim = opts.stim(t);
else
    Istim = opts.stim(:)';
    if numel(Istim) ~= N
        error('simulate_PB_reduced:stimLength', ...
            'Numeric stim must have %d samples to match 0:dt:tEnd (got %d).', N, numel(Istim));
    end
end

V = zeros(1,N); m = zeros(1,N); v = zeros(1,N); f = zeros(1,N); to = zeros(1,N); X = zeros(1,N);
% NaN- not zero-fill when traces are off: see the note in SIMULATE_PB_FULL.
if opts.returnCurrents, blank = @() zeros(1,N); else, blank = @() nan(1,N); end
INa=blank(); ICa=blank(); Ito=blank(); IK=blank(); IK1=blank();
INab=blank(); ICab=blank(); INaK=blank(); INaCa=blank();

V(1) = opts.V0;
g0 = pb_reduced_gates(V(1), params);
m(1)=g0.minf; v(1)=g0.vinf; f(1)=g0.finf; to(1)=g0.toinf; X(1)=g0.Xinf;

if opts.method == "exponential"
    for i = 1:N-1
        g = pb_reduced_gates(V(i), params);
        % Reuse g: pb_reduced_currents would otherwise recompute every
        % gate a second time on every step just to get dinf/rinf.
        c = pb_reduced_currents(V(i), m(i), v(i), f(i), to(i), X(i), params, const, g);

        % Only the STORES are output-only work; the pb_reduced_currents
        % call itself is required by the integration (it supplies EK,
        % ENa, INaK, INaCa below), so it is never skipped.
        if opts.returnCurrents
            INa(i)=c.INa; ICa(i)=c.ICa; Ito(i)=c.Ito; IK(i)=c.IK; IK1(i)=c.IK1;
            INab(i)=c.INab; ICab(i)=c.ICab; INaK(i)=c.INaK; INaCa(i)=c.INaCa;
        end

        gNa_eff = params.gNa*m(i)^3*v(i)^2;
        gCa_eff = params.gCa*g.dinf*f(i)*params.fCa;
        gto_eff = params.gto*g.rinf*to(i);
        gK_eff  = params.gK*X(i)^2;
        aK1 = 0.1/(1+exp(0.06*(V(i)-c.EK-200)));
        bK1 = (3*exp(2e-4*(V(i)-c.EK+100))+exp(0.1*(V(i)-c.EK-10)))/(1+exp(-0.5*(V(i)-c.EK)));
        gK1_eff = params.gK1*(aK1/(aK1+bK1));

        gTot = gNa_eff+gCa_eff+gto_eff+gK_eff+gK1_eff+params.gNab+params.gCab;
        ERevSum = gNa_eff*c.ENa+gCa_eff*c.ECa+gto_eff*c.Eto+gK_eff*c.EK+gK1_eff*c.EK ...
                  +params.gNab*c.ENa+params.gCab*c.ECa;

        Vinf = (ERevSum - c.INaK - c.INaCa + Istim(i)) / gTot;
        tauV = 1/gTot;

        m(i+1)  = g.minf  + (m(i)-g.minf)  *exp(-dt/g.taum);
        v(i+1)  = g.vinf  + (v(i)-g.vinf)  *exp(-dt/g.tauv);
        f(i+1)  = g.finf  + (f(i)-g.finf)  *exp(-dt/g.tauf);
        to(i+1) = g.toinf + (to(i)-g.toinf)*exp(-dt/g.tauto);
        X(i+1)  = g.Xinf  + (X(i)-g.Xinf)  *exp(-dt/g.tauX);
        V(i+1)  = Vinf    + (V(i)-Vinf)    *exp(-dt/tauV);
    end
    if opts.returnCurrents
        % Guarded: writing sample N unconditionally left an array that was
        % entirely unfilled except for its final element.
        c = pb_reduced_currents(V(N), m(N), v(N), f(N), to(N), X(N), params, const);
        INa(N)=c.INa; ICa(N)=c.ICa; Ito(N)=c.Ito; IK(N)=c.IK; IK1(N)=c.IK1;
        INab(N)=c.INab; ICab(N)=c.ICab; INaK(N)=c.INaK; INaCa(N)=c.INaCa;
    end
else
    for i = 1:N-1
        y = [V(i); m(i); v(i); f(i); to(i); X(i)];
        k1 = pb_reduced_derivatives(y,          Istim(i), params, const);
        k2 = pb_reduced_derivatives(y+0.5*dt*k1, Istim(i), params, const);
        k3 = pb_reduced_derivatives(y+0.5*dt*k2, Istim(i), params, const);
        k4 = pb_reduced_derivatives(y+dt*k3,     Istim(i), params, const);
        ynext = y + (dt/6)*(k1+2*k2+2*k3+k4);
        V(i+1)=ynext(1); m(i+1)=ynext(2); v(i+1)=ynext(3);
        f(i+1)=ynext(4); to(i+1)=ynext(5); X(i+1)=ynext(6);
    end
    if opts.returnCurrents
        % The RK4 path previously ignored returnCurrents entirely, so timing
        % this path was never comparable with the exponential one.
        for i = 1:N
            c = pb_reduced_currents(V(i), m(i), v(i), f(i), to(i), X(i), params, const);
            INa(i)=c.INa; ICa(i)=c.ICa; Ito(i)=c.Ito; IK(i)=c.IK; IK1(i)=c.IK1;
            INab(i)=c.INab; ICab(i)=c.ICab; INaK(i)=c.INaK; INaCa(i)=c.INaCa;
        end
    end
end

out = struct('t',t, 'V',V, 'm',m, 'v',v, 'f',f, 'to',to, 'X',X, 'Istim',Istim, ...
             'INa',INa, 'ICa',ICa, 'Ito',Ito, 'IK',IK, 'IK1',IK1, ...
             'INab',INab, 'ICab',ICab, 'INaK',INaK, 'INaCa',INaCa, ...
             'params',params, 'dt',dt, 'method',char(opts.method), 'model','PB_reduced');
end
