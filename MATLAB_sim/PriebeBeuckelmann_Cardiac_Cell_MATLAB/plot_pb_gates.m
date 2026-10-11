function fig = plot_pb_gates(res)
%   For the reduced model the two adiabatically eliminated gates are
%   overlaid as dotted curves: d and r were fast enough that Bernus et al.
%   replaced them by their instantaneous steady states d_inf(V) and
%   r_inf(V) (PB_REDUCED_GATES), so they carry no state of their own, and
%   the single gate v (entering as v^2) replaces the h*j inactivation
%   product of the full model. Seeing d_inf(V) tracking the voltage
%   instantaneously, next to the genuinely slow f, is the clearest single
%   picture of what the reduction actually did.
%   Every gate is dimensionless and bounded to [0 1], so all panels share
%   ylim([0 1]).
arguments
    res (1,1) struct
end

pal = pb_palette();
isFull = localIsFull(res);

[k, tWin, ap] = localPickBeat(res);
idx = res.t >= tWin(1) & res.t <= tWin(2);
if k > 0
    t0 = ap.stimTimes(k);
else
    t0 = 0;
end
tt = res.t(idx) - t0;
Vseg = res.V(idx);
Vnorm = (Vseg - min(Vseg)) / max(max(Vseg) - min(Vseg), eps);

white = [1 1 1];
colJ   = 0.60*pal.v   + 0.40*white;
colR   = 0.55*pal.Ito + 0.45*white;
colXr  = 0.60*pal.X   + 0.40*white;
colFCa = 0.55*pal.f   + 0.45*white;

fig = figure('Color','w','Position',[100 100 960 820]);
tl = tiledlayout(4, 1, 'TileSpacing','compact');

if isFull
    localTile(tt, Vnorm, pal);
    plot(tt, res.m(idx), '-',  'Color', pal.m, 'LineWidth', 1.4, 'DisplayName','m (activation)');
    plot(tt, res.h(idx), '-',  'Color', pal.v, 'LineWidth', 1.4, 'DisplayName','h (fast inactivation)');
    plot(tt, res.j(idx), '--', 'Color', colJ,  'LineWidth', 1.4, 'DisplayName','j (slow inactivation)');
    localFinishTile('I_{Na} gates (current \propto m^3 h j)');

    localTile(tt, Vnorm, pal);
    plot(tt, res.d(idx), '-', 'Color', pal.d, 'LineWidth', 1.4, 'DisplayName','d (activation)');
    plot(tt, res.f(idx), '-', 'Color', pal.f, 'LineWidth', 1.4, 'DisplayName','f (voltage inactivation)');
    % f_Ca is NOT a state variable: it is the algebraic Ca2+-dependent
    % inactivation factor 1/(1+[Ca2+]_i/K_m) that multiplies d*f in the
    % full model. 
    if isfield(res,'Cai') && isfield(res,'params') && isfield(res.params,'KmCa_fCa')
        fCa = 1 ./ (1 + res.Cai(idx)/res.params.KmCa_fCa);
        plot(tt, fCa, ':', 'Color', colFCa, 'LineWidth', 1.4, ...
             'DisplayName','f_{Ca}([Ca^{2+}]_i) (algebraic, not a gate)');
    end
    localFinishTile('I_{Ca} gates (current \propto d f f_{Ca})');

    localTile(tt, Vnorm, pal);
    plot(tt, res.r(idx),  '--', 'Color', colR,   'LineWidth', 1.4, 'DisplayName','r (activation)');
    plot(tt, res.to(idx), '-',  'Color', pal.to, 'LineWidth', 1.4, 'DisplayName','t_o (inactivation)');
    localFinishTile('I_{to} gates (current \propto r t_o)');

    localTile(tt, Vnorm, pal);
    plot(tt, res.Xs(idx), '-',  'Color', pal.X, 'LineWidth', 1.4, 'DisplayName','X_s (I_{Ks} activation)');
    plot(tt, res.Xr(idx), '--', 'Color', colXr, 'LineWidth', 1.4, 'DisplayName','X_r (I_{Kr} activation)');
    localFinishTile('Delayed-rectifier gates (I_{Ks} \propto X_s^2, I_{Kr} \propto X_r)');

    noteLines = {};
else
    haveElim = isfield(res,'params') && all(isfield(res.params, {'Vshift','invP','cellType'}));
    if haveElim
        g = pb_reduced_gates(Vseg, res.params);
    end

    localTile(tt, Vnorm, pal);
    plot(tt, res.m(idx), '-', 'Color', pal.m, 'LineWidth', 1.4, 'DisplayName','m (activation)');
    plot(tt, res.v(idx), '-', 'Color', pal.v, 'LineWidth', 1.4, ...
         'DisplayName','v (single inactivation gate, replaces h j)');
    localFinishTile('I_{Na} gates (current \propto m^3 v^2)');

    localTile(tt, Vnorm, pal);
    if haveElim
        plot(tt, g.dinf, ':', 'Color', pal.d, 'LineWidth', 1.4, ...
             'DisplayName','d_{\infty}(V) (eliminated: instantaneous, no state)');
    end
    plot(tt, res.f(idx), '-', 'Color', pal.f, 'LineWidth', 1.4, 'DisplayName','f (voltage inactivation)');
    if isfield(res,'params') && isfield(res.params,'fCa')
        plot(tt, res.params.fCa*ones(size(tt)), ':', 'Color', colFCa, 'LineWidth', 1.4, ...
             'DisplayName', sprintf('f_{Ca} = %.2f (frozen: [Ca^{2+}]_i constant)', res.params.fCa));
    end
    localFinishTile('I_{Ca} gates (current \propto d_{\infty}(V) f f_{Ca})');

    localTile(tt, Vnorm, pal);
    if haveElim
        plot(tt, g.rinf, ':', 'Color', colR, 'LineWidth', 1.4, ...
             'DisplayName','r_{\infty}(V) (eliminated: instantaneous, no state)');
    end
    plot(tt, res.to(idx), '-', 'Color', pal.to, 'LineWidth', 1.4, 'DisplayName','t_o (inactivation)');
    localFinishTile('I_{to} gates (current \propto r_{\infty}(V) t_o)');

    localTile(tt, Vnorm, pal);
    plot(tt, res.X(idx), '-', 'Color', pal.X, 'LineWidth', 1.4, ...
         'DisplayName','X (single merged delayed rectifier)');
    localFinishTile('Delayed-rectifier gate (I_K \propto X^2, one gate for I_{Kr}+I_{Ks})');

    noteLines = { ...
        'Reduction: d and r are adiabatically eliminated - they survive only as the instantaneous steady states d_{\infty}(V), r_{\infty}(V) (dotted),', ...
        'and the single gate v (entering as v^2) replaces the full model''s h j inactivation product: 10 gates become 5.'};
end

xlabel('Time after stimulus (ms)');

if k > 0
    sub = sprintf('gating variables, beat %d of %d, aligned on its stimulus (t = 0)', k, ap.count);
else
    sub = 'gating variables, whole trace shown: no paced beats detected';
end
title(tl, [{localModelLabel(res)}, {sub}, noteLines]);
end


function localTile(tt, Vnorm, pal)
nexttile;
plot(tt, Vnorm, ':', 'Color', pal.stim, 'LineWidth', 1.0, ...
     'DisplayName','V_m (normalized, reference only)');
hold on;
end


function localFinishTile(ttl)
ylim([0 1]);
ylabel('Gate value');
title(ttl);
grid on; box off;
legend('Location','best', 'Box','off', 'FontSize', 8);
end


function tf = localIsFull(res)
if ~isfield(res, 'model')
    error('plot_pb_gates:noModelField', ...
        'RES has no "model" field; expected the output of simulate_PB_full or simulate_PB_reduced.');
end
switch res.model
    case 'PB_full',    tf = true;
    case 'PB_reduced', tf = false;
    otherwise
        error('plot_pb_gates:unknownModel', ...
            'Unrecognized res.model "%s" (expected "PB_full" or "PB_reduced").', res.model);
end
end


function [k, tWin, ap] = localPickBeat(res)
ap = detect_pb_ap(res);
if ap.count == 0
    k = 0;
    tWin = [res.t(1) res.t(end)];
    return
end
good = find(isfinite(ap.APD(:,end)));
if isempty(good), k = ap.count; else, k = good(end); end
if k < ap.count
    tStop = ap.stimTimes(k+1);     
else
    tStop = res.t(end);
end
pre  = min(20, ap.stimTimes(k) - res.t(1));  
tWin = [ap.stimTimes(k)-pre, tStop];
end


function s = localModelLabel(res)
if isfield(res,'model'), mdl = res.model; else, mdl = ''; end
switch mdl
    case 'PB_full'
        s = 'Priebe-Beuckelmann (1998) full model, 15 state variables';
    case 'PB_reduced'
        s = 'Bernus et al. (2002) reduced model, 6 state variables';
    otherwise
        s = 'Priebe-Beuckelmann human ventricular myocyte';
end
if isfield(res,'params')
    if isfield(res.params,'condition')
        s = sprintf('%s - %s myocyte', s, char(res.params.condition));
    elseif isfield(res.params,'cellType')
        s = sprintf('%s - %s cell', s, char(res.params.cellType));
    end
end
end
