function fig = plot_pb_currents(res)
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

fig = figure('Color','w','Position',[100 100 1000 940]);
tl = tiledlayout(5, 2, 'TileSpacing','compact');

axV = nexttile([1 2]); hold on;
plot(tt, res.V(idx), '-', 'Color', pal.V, 'LineWidth', 1.6);
localFinishPanel(pal, 'Membrane potential', 'V_m (mV)', false, false);

axNa = nexttile; hold on;
plot(tt, res.INa(idx), '-', 'Color', pal.INa, 'LineWidth', 1.3);
localFinishPanel(pal, 'I_{Na}: whole beat', 'I_{Na} (pA/pF)', true, false);

axNaZoom = nexttile; hold on;
plot(tt, res.INa(idx), '-', 'Color', pal.INa, 'LineWidth', 1.3);
localFinishPanel(pal, 'I_{Na}: upstroke detail', 'I_{Na} (pA/pF)', true, false);
xlim(axNaZoom, [max(tt(1),-2) min(tt(end),6)]);
peakINa = min(res.INa(idx));        
txt = text(0, 0, sprintf('peak %.0f pA/pF', peakINa), 'Parent', axNaZoom, ...
           'HorizontalAlignment','right', 'FontSize', 8.5);
txt.Units = 'normalized';           
txt.Position = [0.96 0.12 0];        
axCa = nexttile([1 2]); hold on;
plot(tt, res.ICa(idx), '-', 'Color', pal.ICa, 'LineWidth', 1.3, 'DisplayName','I_{Ca} (L-type)');
plot(tt, res.Ito(idx), '-', 'Color', pal.Ito, 'LineWidth', 1.3, 'DisplayName','I_{to}');
localFinishPanel(pal, 'Plateau-shaping currents', 'Current (pA/pF)', true, true);

axK = nexttile([1 2]); hold on;
if isFull
    colIKs = 0.60*pal.IK + 0.40*[1 1 1];
    colSum = 0.60*pal.IK;
    plot(tt, res.IKr(idx), '-',  'Color', pal.IK,  'LineWidth', 1.3, 'DisplayName','I_{Kr}');
    plot(tt, res.IKs(idx), '--', 'Color', colIKs,  'LineWidth', 1.3, 'DisplayName','I_{Ks}');
    plot(tt, res.IKr(idx)+res.IKs(idx), ':', 'Color', colSum, 'LineWidth', 1.1, ...
         'DisplayName','I_{Kr}+I_{Ks} (= reduced model''s I_K)');
    plot(tt, res.IK1(idx), '-', 'Color', pal.IK1, 'LineWidth', 1.3, 'DisplayName','I_{K1}');
    ttlK = 'K^+ currents: delayed rectifier resolved into I_{Kr} and I_{Ks}';
else
    plot(tt, res.IK(idx),  '-', 'Color', pal.IK,  'LineWidth', 1.3, 'DisplayName','I_K (merged I_{Kr}+I_{Ks})');
    plot(tt, res.IK1(idx), '-', 'Color', pal.IK1, 'LineWidth', 1.3, 'DisplayName','I_{K1}');
    ttlK = 'K^+ currents: single merged delayed rectifier I_K';
end
localFinishPanel(pal, ttlK, 'Current (pA/pF)', true, true);

axP = nexttile([1 2]); hold on;
plot(tt, res.INaK(idx),  '-', 'Color', pal.INaK,  'LineWidth', 1.3, 'DisplayName','I_{NaK} (Na^+/K^+ pump)');
plot(tt, res.INaCa(idx), '-', 'Color', pal.INaCa, 'LineWidth', 1.3, 'DisplayName','I_{NaCa} (Na^+/Ca^{2+} exchanger)');
localFinishPanel(pal, 'Electrogenic transporters', 'Current (pA/pF)', true, true);
xlabel('Time after stimulus (ms)');

linkaxes([axV axNa axCa axK axP], 'x');
xlim(axV, [tt(1) tt(end)]);

if k > 0
    sub = sprintf('beat %d of %d, aligned on its stimulus (t = 0)', k, ap.count);
else
    sub = 'whole trace shown: no paced beats detected';
end
title(tl, {localModelLabel(res), ['ionic currents, ' sub]});
end


function localFinishPanel(pal, ttl, ylab, zeroLine, showLegend)
if zeroLine
    yline(0, '-', 'Color', pal.grid, 'LineWidth', 0.8, 'HandleVisibility','off');
end
xline(0, ':', 'Color', pal.stim, 'LineWidth', 0.5, 'HandleVisibility','off');
title(ttl);
ylabel(ylab);
grid on; box off;
if showLegend
    legend('Location','best', 'Box','off');
end
end


function tf = localIsFull(res)
if ~isfield(res, 'model')
    error('plot_pb_currents:noModelField', ...
        'RES has no "model" field; expected the output of simulate_PB_full or simulate_PB_reduced.');
end
switch res.model
    case 'PB_full',    tf = true;
    case 'PB_reduced', tf = false;
    otherwise
        error('plot_pb_currents:unknownModel', ...
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
