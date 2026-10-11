function worst = verify_pb_rhs(throwOnFail)
%VERIFY_PB_RHS Guard the optimization in PB_FULL_RHS.
%   PB_FULL_RHS must be numerically identical to the readable reference
%   PB_FULL_DERIVATIVES, because only the latter is transcribed
%   one-to-one from the Priebe-Beuckelmann CellML source; the former is a
%   hand-inlined, allocation-free rewrite kept for speed. Any divergence
%   means the fast path has drifted from the equations of record.
%
%   WORST = VERIFY_PB_RHS() compares the two across a spread of
%   physiologically reachable states (rest, upstroke, plateau with CICR
%   active, repolarization, plus 8 randomized states), for both stimulus
%   levels and BOTH parameter conditions, and returns the largest
%   relative difference found.
%
%   VERIFY_PB_RHS(TRUE) additionally raises an error when the tolerance
%   is exceeded, so an automated caller (SMOKE_TEST) can actually fail on
%   a regression. A check that only prints "FAIL" guards nothing.

if nargin < 1, throwOnFail = false; end
worst = 0;
for cond = ["nonfailing", "failing"]
    worst = max(worst, local_check(cond, throwOnFail));
end
end

% -------------------------------------------------------------------------
function worst = local_check(condition, throwOnFail)
p = pb_full_params(condition);

y0 = [p.V0;p.m0;p.h0;p.j0;p.d0;p.f0;p.r0;p.to0;p.Xs0;p.Xr0; ...
      p.Nai0;p.Ki0;p.Cai0;p.CaJSR0;p.CaNSR0;0;0;0];

states = {};
states{end+1} = y0;                                              %#ok<*AGROW> rest
yv = y0; yv(1) = -35; yv(2) = 0.6; yv(16) = 0.9;  states{end+1} = yv;   % upstroke
yv = y0; yv(1) = 10;  yv(5) = 0.5; yv(6) = 0.5; yv(13) = 0.0007; ...
         yv(14) = 1.2; yv(16) = 0.15; yv(17) = 0.5; yv(18) = 1e-4;
states{end+1} = yv;                                              % plateau + CICR active
yv = y0; yv(1) = -60; yv(9) = 0.3; yv(10) = 0.4; yv(13) = 0.0004;
states{end+1} = yv;                                              % repolarizing

rng(0);
for k = 1:8
    yv = y0;
    yv(1) = -95 + 140*rand;                 % V in [-95, 45]
    yv(2:10) = rand(9,1);                   % gates in [0,1]
    yv(11) = 8 + 4*rand;                    % Nai
    yv(12) = 135 + 10*rand;                 % Ki
    yv(13) = 1e-4 + 1e-3*rand;              % Cai
    yv(14) = 0.5 + 2*rand;                  % CaJSR
    yv(15) = 0.5 + 2*rand;                  % CaNSR
    yv(16) = rand; yv(17) = rand; yv(18) = 2e-4*rand;
    states{end+1} = yv;
end

stims  = [0, 15];
maxRel = zeros(18,1);
worst  = 0;
for s = 1:numel(states)
    for q = 1:numel(stims)
        a = pb_full_derivatives(states{s}, stims(q), p);
        b = pb_full_rhs(states{s},         stims(q), p);
        den = max(abs(a), 1e-12);
        rel = abs(a-b)./den;
        maxRel = max(maxRel, rel);
        worst  = max(worst, max(rel));
    end
end

% Tolerance note: the two implementations are ALGEBRAICALLY identical but
% not bit-identical, because the reference evaluates each gate as
% (xinf - x)/tau with xinf = a/(a+b), tau = 1/(a+b), whereas the fast
% path evaluates the expanded-and-cancelled form a*(1-x) - b*x. Those
% agree exactly in exact arithmetic and differ only by floating-point
% rounding (~1e-11 relative on the gates with the largest a+b). A
% threshold of 1e-12 would flag ordinary round-off as a bug; 1e-8 is
% still ~7 orders of magnitude below any physically meaningful
% difference, while remaining tight enough to catch a genuine
% transcription error (which would show up as O(1e-2) or larger).
tol = 1e-8;
fprintf('   [%-10s] worst relative difference %.3e  (tol %.0e)\n', condition, worst, tol);
if worst >= tol
    msg = sprintf(['PB_FULL_RHS diverges from PB_FULL_DERIVATIVES for the "%s" ' ...
                   'parameter set (worst relative difference %.3e > tol %.0e). ' ...
                   'The optimized path has a transcription bug.'], condition, worst, tol);
    fprintf('   FAIL - %s\n', msg);
    if throwOnFail
        error('verify_pb_rhs:mismatch', '%s', msg);
    end
end
end
