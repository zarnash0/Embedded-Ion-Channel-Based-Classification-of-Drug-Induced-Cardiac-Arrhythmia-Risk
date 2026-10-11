function out = simulate_PB_full(params, tEnd, opts)
%SIMULATE_PB_FULL Integrate the full (15-variable, "second-generation")
%   Priebe-Beuckelmann (1998) human ventricular cell model.
%   OUT = SIMULATE_PB_FULL(PARAMS, TEND) integrates from t=0 to t=TEND
%   (ms), paced by a periodic stimulus train (default: 1 Hz, 3-ms pulses
%   of 15 pA/pF, matching the CellML reference protocol -- see
%   SIMULATE_PB_REDUCED for the sign-convention note).
%
%   Name-value options: dt (default 0.01 ms), stim, V0 (mV; the other 17
%   initial state values are always taken from PARAMS, matching the
%   CellML reference exactly).
%
%   INTEGRATION: classical 4th-order Runge-Kutta with a fixed time step,
%   exactly the method Priebe & Beuckelmann themselves report ("A
%   fourth-order Runge-Kutta method with fixed time intervals is used for
%   numerical integration of differential equations"). Unlike the
%   Hodgkin-Huxley/Noble/reduced-PB models earlier in this project, this
%   full model's ion-concentration and Ca2+-handling state variables are
%   not simple relaxation (Rush-Larsen-amenable) equations, so RK4 on the
%   raw system -- rather than an exponential/semi-implicit scheme -- is
%   both the simplest and the most source-faithful choice here.

arguments
    params (1,1) struct
    tEnd (1,1) double {mustBePositive}
    opts.dt (1,1) double {mustBePositive} = 0.005
    opts.stim = []
    opts.V0 (1,1) double = params.V0
    opts.returnCurrents (1,1) logical = true
end
% returnCurrents=false skips the post-integration pass that rebuilds the
% individual current traces. That pass is O(N) and costs a large,
% model-specific fraction of the total, so leaving it inside a timed
% region would make a full-vs-reduced speed comparison measure
% bookkeeping rather than integration. COMPARE_PB_MODELS therefore times
% both models with returnCurrents=false. Physics is unaffected.

t = 0:opts.dt:tEnd;
N = numel(t);
dt = opts.dt;

if isempty(opts.stim)
    % 50 ms offset so each beat is preceded by a visible resting
    % baseline (needed for a clean RMP measurement, and it gives every
    % stimulus a detectable rising edge).
    stimFcn = pb_make_train(1000, 3, 15, 50);
    Istim = stimFcn(t);
elseif isa(opts.stim, 'function_handle')
    Istim = opts.stim(t);
else
    Istim = opts.stim(:)';
    if numel(Istim) ~= N
        error('simulate_PB_full:stimLength', ...
            'Numeric stim must have %d samples to match 0:dt:tEnd (got %d).', N, numel(Istim));
    end
end

y = zeros(18, N);
y(:,1) = [opts.V0; params.m0; params.h0; params.j0; params.d0; params.f0; params.r0; params.to0; ...
          params.Xs0; params.Xr0; params.Nai0; params.Ki0; params.Cai0; params.CaJSR0; params.CaNSR0; ...
          params.APtrack0; params.APtrack20; params.Cainfluxtrack0];

% PB_FULL_RHS is the allocation-free twin of PB_FULL_DERIVATIVES (the two
% are checked against each other by VERIFY_PB_RHS); using it here is a
% pure speed measure, ~10x faster because it avoids per-substage struct
% construction.
for i = 1:N-1
    yi = y(:,i);
    k1 = pb_full_rhs(yi,            Istim(i), params);
    k2 = pb_full_rhs(yi+0.5*dt*k1,  Istim(i), params);
    k3 = pb_full_rhs(yi+0.5*dt*k2,  Istim(i), params);
    k4 = pb_full_rhs(yi+dt*k3,      Istim(i), params);
    y(:,i+1) = yi + (dt/6)*(k1+2*k2+2*k3+k4);
end

% Zero-fill only when these are actually going to be filled. With
% returnCurrents=false the arrays are returned under their normal names but
% never written, and an all-zero current trace looks like legitimate data --
% a caller who forgot the flag would plot and analyse it without noticing.
% NaN makes that failure loud. (Same reasoning applied in SIMULATE_HH and
% SIMULATE_NOBLE.)
if opts.returnCurrents, blank = @() zeros(1,N); else, blank = @() nan(1,N); end
iNa=blank(); iCa=blank(); ito=blank(); iKr=blank(); iKs=blank(); iK1=blank();
iNab=blank(); iCab=blank(); iNaK=blank(); iNaCa=blank();
iup=blank(); ileak=blank(); itr=blank();
if opts.returnCurrents
    for i = 1:N
        c = pb_full_currents(y(1:15,i), params);
        iNa(i)=c.iNa; iCa(i)=c.iCa; ito(i)=c.ito; iKr(i)=c.iKr; iKs(i)=c.iKs; iK1(i)=c.iK1;
        iNab(i)=c.iNab; iCab(i)=c.iCab; iNaK(i)=c.iNaK; iNaCa(i)=c.iNaCa;
        iup(i)=c.iup; ileak(i)=c.ileak; itr(i)=c.itr;
    end
end

out = struct('t',t, 'V',y(1,:), 'm',y(2,:), 'h',y(3,:), 'j',y(4,:), 'd',y(5,:), 'f',y(6,:), ...
    'r',y(7,:), 'to',y(8,:), 'Xs',y(9,:), 'Xr',y(10,:), ...
    'Nai',y(11,:), 'Ki',y(12,:), 'Cai',y(13,:), 'CaJSR',y(14,:), 'CaNSR',y(15,:), ...
    'APtrack',y(16,:), 'APtrack2',y(17,:), 'Cainfluxtrack',y(18,:), ...
    'Istim',Istim, 'INa',iNa, 'ICa',iCa, 'Ito',ito, 'IKr',iKr, 'IKs',iKs, 'IK1',iK1, ...
    'INab',iNab, 'ICab',iCab, 'INaK',iNaK, 'INaCa',iNaCa, 'Iup',iup, 'Ileak',ileak, 'Itr',itr, ...
    'params',params, 'dt',dt, 'method','rk4', 'model','PB_full');
end
