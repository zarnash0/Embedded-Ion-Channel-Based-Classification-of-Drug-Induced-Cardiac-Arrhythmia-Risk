function fig = plot_pb_comparison(cmp)
%     A  One representative beat, both models overlaid -- the reduction
%        preserves AP shape (mirrors Bernus et al. 2002, Fig. 3A).
%     B  Beat-to-beat APD90 for both models over the whole run -- the full
%        model's APD creeps while the reduced model's is flat (mirrors
%        Bernus Figs. 9A / 10C).
%     C  The full model's [Na+]i and [K+]i over the run, against the fixed
%        10 / 140 mM the reduced model uses by construction (mirrors
%        Bernus Figs. 9B, 9C).
%     D  The efficiency numbers: wall-clock time and state-variable count.
arguments
    cmp (1,1) struct
end

pal = pb_palette();
p = cmp.protocol;

labFull = sprintf('full, 15 vars (dt = %g ms)', cmp.full.dt);
labRed  = sprintf('reduced, 6 vars (dt = %g ms)', cmp.reduced.dt);

fig = figure('Color','w','Position',[70 55 1220 850]);
tl = tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');
ax = nexttile;
hold(ax, 'on');
preMs = 20;
tStimF = cmp.apFull.stimTimes(cmp.beatFull);
tStimR = cmp.apReduced.stimTimes(cmp.beatReduced);
[tF, vF] = localBeatWindow(cmp.full,    tStimF, preMs, p.searchWindow);
[tR, vR] = localBeatWindow(cmp.reduced, tStimR, preMs, p.searchWindow);
hF = plot(ax, tF, vF, 'Color', pal.full,    'LineWidth', 1.8, 'DisplayName', labFull);
hR = plot(ax, tR, vR, 'Color', pal.reduced, 'LineWidth', 1.8, 'DisplayName', labRed);
yline(0, ':', 'Color', pal.grid, 'HandleVisibility','off');   
xlabel(ax, 'Time after stimulus (ms)');
ylabel(ax, 'V_m (mV)');
title(ax, sprintf('A   AP shape survives the reduction (beat %d)', cmp.beatFull), ...
    'FontWeight','bold');
legend(ax, [hF hR], 'Location','northeast', 'Box','off');
xlim(ax, [-preMs p.searchWindow]);
localPadYLim(ax, 0.24, 0.08);
localTidy(ax, pal);

a90  = localMetric(cmp, 'APD90');
adv  = localMetric(cmp, 'dV/dt max');
armp = localMetric(cmp, 'RMP');
txt = { sprintf('%-16s%8s%9s%9s', '', 'full', 'reduced', '|diff|'), ...
        sprintf('%-16s%8.1f%9.1f%9.1f', 'APD90 (ms)',      a90(1),  a90(2),  a90(3)), ...
        sprintf('%-16s%8.1f%9.1f%9.1f', 'dV/dt max (V/s)', adv(1),  adv(2),  adv(3)), ...
        sprintf('%-16s%8.1f%9.1f%9.1f', 'RMP (mV)',        armp(1), armp(2), armp(3)) };
text(ax, 0.02, 0.03, txt, 'Units','normalized', 'FontSize', 9, ...
    'VerticalAlignment','bottom', 'FontName','FixedWidth', 'Color', pal.V);
ax = nexttile;
hold(ax, 'on');
hF = plot(ax, cmp.drift.beatsFull, cmp.drift.APD90Full, '-o', ...
    'Color', pal.full, 'MarkerFaceColor', pal.full, 'MarkerSize', 5, ...
    'LineWidth', 1.6, 'DisplayName', labFull);
hR = plot(ax, cmp.drift.beatsReduced, cmp.drift.APD90Reduced, '-s', ...
    'Color', pal.reduced, 'MarkerFaceColor', pal.reduced, 'MarkerSize', 5, ...
    'LineWidth', 1.6, 'DisplayName', labRed);
xlabel(ax, 'Beat number');
ylabel(ax, 'APD_{90} (ms)');
title(ax, 'B   Beat-to-beat stability', 'FontWeight','bold');
subtitle(ax, sprintf('drift across the run:  full %+.1f ms,   reduced %+.1f ms', ...
    cmp.drift.APD90ChangeFull, cmp.drift.APD90ChangeReduced), ...
    'FontSize', 9, 'FontWeight','normal');
legend(ax, [hF hR], 'Location','northeast', 'Box','off');
nB = max([numel(cmp.drift.beatsFull), numel(cmp.drift.beatsReduced), 2]);
xlim(ax, [0.5, nB+0.5]);
xticks(ax, 1:nB);
localPadYLim(ax, 0.32, 0.30);
localTidy(ax, pal);
text(ax, 0.02, 0.03, sprintf(['Bernus et al.: the full model paced 1000 s at 1 Hz never reaches\n' ...
    'steady state (APD_{90} creeps ~360 \\rightarrow 425 ms as [Na^{+}]_i falls);\n' ...
    'the reduced model has no slow variable to creep with.']), ...
    'Units','normalized', 'FontSize', 8.5, 'VerticalAlignment','bottom', 'Color', pal.V);
ax = nexttile;
stride = max(1, floor(numel(cmp.full.t)/20000));
tC = cmp.full.t(1:stride:end)/1000;  

yyaxis(ax, 'left');
hNa = plot(ax, tC, cmp.full.Nai(1:stride:end), '-', 'Color', pal.INa, 'LineWidth', 1.8, ...
    'DisplayName', '[Na^{+}]_i, full model');
hold(ax, 'on');
hFix = yline(cmp.drift.NaiFixedReduced, '--', 'Color', pal.reduced, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('reduced model: fixed at %g / %g mM', ...
    cmp.drift.NaiFixedReduced, cmp.drift.KiFixedReduced));
ylabel(ax, '[Na^{+}]_i (mM)');
ax.YAxis(1).Color = pal.INa;
localPadYLim(ax, 0.10, 0.45);   

yyaxis(ax, 'right');
hK = plot(ax, tC, cmp.full.Ki(1:stride:end), '-', 'Color', pal.IK, 'LineWidth', 1.8, ...
    'DisplayName', '[K^{+}]_i, full model');
hold(ax, 'on');
yline(cmp.drift.KiFixedReduced, '--', 'Color', pal.reduced, 'LineWidth', 1.5, ...
    'HandleVisibility','off');
ylabel(ax, '[K^{+}]_i (mM)');
ax.YAxis(2).Color = pal.IK;
localPadYLim(ax, 0.10, 0.45);

xlabel(ax, 'Time (s)');
title(ax, 'C   Only the full model has anywhere to drift to', 'FontWeight','bold');
subtitle(ax, sprintf(['over %g s: [Na^{+}]_i %+.4f mM, [K^{+}]_i %+.4f mM ' ...
    '(note the expanded scale)'], cmp.opts.tEnd/1000, cmp.drift.NaiChange, cmp.drift.KiChange), ...
    'FontSize', 9, 'FontWeight','normal');
legend(ax, [hNa hK hFix], 'Location','northwest', 'Box','off');
localTidy(ax, pal);
ax = nexttile;
ratios = [1, cmp.timeReduced/cmp.timeFull; ...
          1, cmp.timeReducedSameDt/cmp.timeFull; ...
          1, cmp.nStateReduced/cmp.nStateFull];
absFull = { sprintf('%.1f s', cmp.timeFull), sprintf('%.1f s', cmp.timeFull), ...
            sprintf('%d', cmp.nStateFull) };
absRed  = { sprintf('%.1f s', cmp.timeReduced), sprintf('%.1f s', cmp.timeReducedSameDt), ...
            sprintf('%d', cmp.nStateReduced) };

b = bar(ax, ratios, 'grouped', 'EdgeColor','none');
b(1).FaceColor = pal.full;     b(1).DisplayName = labFull;
b(2).FaceColor = pal.reduced;  b(2).DisplayName = labRed;
hold(ax, 'on');
text(ax, b(1).XEndPoints, b(1).YEndPoints, absFull, 'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', 'FontSize', 9, 'Color', pal.full);
text(ax, b(2).XEndPoints, b(2).YEndPoints, absRed, 'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', 'FontSize', 9, 'Color', pal.reduced);

xticks(ax, 1:3);
xticklabels(ax, { sprintf('wall clock\neach at its own dt'), ...
                  sprintf('wall clock\nboth at dt = %g ms', cmp.full.dt), ...
                  sprintf('state variables\nper cell') });
ylabel(ax, 'Cost relative to the full model');
ylim(ax, [0 1.32]);
title(ax, 'D   What the reduction costs to run', 'FontWeight','bold');
subtitle(ax, sprintf('%.1fx faster in practice,  %.1fx at equal dt,  %.1fx fewer state variables', ...
    cmp.speedup, cmp.speedupSameDt, cmp.nStateFull/cmp.nStateReduced), ...
    'FontSize', 9, 'FontWeight','normal');
legend(ax, [b(1) b(2)], 'Location','northeast', 'Box','off');
localTidy(ax, pal);

title(tl, sprintf(['Priebe-Beuckelmann (1998) full model  vs  Bernus et al. (2002) reduced model' ...
    '        %s myocyte / %s cell,  %g Hz pacing,  %g s run'], ...
    cmp.opts.condition, cmp.opts.cellType, 1000/p.period, cmp.opts.tEnd/1000), ...
    'FontWeight','bold', 'FontSize', 12);
end
function [tRel, V] = localBeatWindow(res, tStim, preMs, postMs)
mask = res.t >= (tStim - preMs) & res.t <= (tStim + postMs);
tRel = res.t(mask) - tStim;
V = res.V(mask);
end
function v = localMetric(cmp, name)
row = cmp.apdTable.Metric == string(name);
v = [cmp.apdTable.Full(row), cmp.apdTable.Reduced(row), cmp.apdTable.AbsDiff(row)];
end

function localPadYLim(ax, fracBelow, fracAbove)
yl = ylim(ax);
r = diff(yl);
if r <= 0
    return
end
ylim(ax, [yl(1)-fracBelow*r, yl(2)+fracAbove*r]);
end
function localTidy(ax, pal)
grid(ax, 'on');
box(ax, 'off');
ax.GridColor = pal.grid;
ax.GridAlpha = 1;
ax.FontSize = 10;
ax.Layer = 'top';
end
