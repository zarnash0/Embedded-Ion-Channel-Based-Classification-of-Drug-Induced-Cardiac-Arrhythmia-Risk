function fig = plot_annotated_ap(res, opts)
arguments
    res (1,1) struct
    opts.threshold (1,1) double = -10
end

pal = hh_palette();
sp  = detect_spikes(res.t, res.V, threshold=opts.threshold);

fig = figure('Color','w','Position',[100 100 850 520]);
plot(res.t, res.V, 'Color', pal.V, 'LineWidth', 1.8); hold on;
yline(opts.threshold, ':', 'spike threshold', 'Color', [0.55 0.55 0.55]);
xlabel('Time (ms)'); ylabel('Membrane potential V_m (mV)');
title('Hodgkin-Huxley action potential — annotated phases');
grid on; box off;

if sp.count >= 1
    tCross = sp.times(1);
    tPeak  = sp.peakTimes(1);
    pkIdx  = find(res.t >= tPeak, 1, 'first');
    pkV    = res.V(pkIdx);

    winEnd = min(numel(res.t), pkIdx + round(15/res.dt));
    [undV, relIdx] = min(res.V(pkIdx:winEnd));
    undIdx = pkIdx + relIdx - 1;
    tUnd   = res.t(undIdx);
    tReturn = min(res.t(end), tUnd + 5);

    labels = {'1','2','3','4','5','6'};
    ptsT = [tCross, (tCross+tPeak)/2, tPeak, (tPeak+tUnd)/2, tUnd, tReturn];
    ptsV = interp1(res.t, res.V, ptsT);

    plot(ptsT, ptsV, 'o', 'MarkerFaceColor', [0.85 0.10 0.10], ...
         'MarkerEdgeColor','k', 'MarkerSize', 6, 'LineStyle','none');
    for i = 1:6
        text(ptsT(i), ptsV(i)+4, labels{i}, 'FontWeight','bold', 'FontSize', 11, ...
             'HorizontalAlignment','center');
    end

    legendText = { ...
        '1: threshold reached — Na^+ channels open', ...
        '2: depolarizing upstroke — Na^+ influx', ...
        '3: peak — Na^+ channels inactivate', ...
        '4: repolarization — K^+ efflux dominates', ...
        '5: undershoot — delayed K^+ closure', ...
        '6: return to resting potential'};
    ax = gca;
    text(ax.XLim(2), ax.YLim(1), '', 'HandleVisibility','off');
    annotation('textbox', [0.58 0.15 0.4 0.35], 'String', legendText, ...
        'EdgeColor','none', 'FontSize', 8.5, 'BackgroundColor','w');
end
end
