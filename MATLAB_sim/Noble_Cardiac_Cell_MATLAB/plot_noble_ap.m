function fig = plot_noble_ap(res, opts)
%   Annotates the 2nd detected beat (skipping the very first, which can
%   carry initial-condition transients) with the standard cardiac
%   action-potential phase numbering. This 1962 model predates the
%   transient outward current (I_to) that produces the phase-1 notch in
%   later human-ventricular models
arguments
    res (1,1) struct
    opts.threshold (1,1) double = -40
end

pal = noble_palette();
fig = figure('Color','w','Position',[100 100 950 550]);
plot(res.t, res.V, 'Color', pal.V, 'LineWidth', 1.6); hold on;
yline(opts.threshold, ':', 'Color', [0.6 0.6 0.6]);
xlabel('Time (ms)'); ylabel('Membrane potential V_m (mV)');
title('Noble (1962) spontaneous Purkinje-fibre pacemaker + action potentials');
grid on; box off;
beats = detect_beats(res, threshold=opts.threshold);
if beats.count >= 1
    lastGoodIdx = find(~isnan(beats.APD), 1, 'last');
    if ~isempty(lastGoodIdx) && lastGoodIdx < beats.count
        xlim([0, min(res.t(end), beats.times(lastGoodIdx+1) + 50)]);
    end
end

if beats.count >= 2
    k = 2;
    tCross = beats.times(k);
    tPeak  = beats.peakTimes(k);
    apdEff = beats.APD(k);
    if isnan(apdEff), apdEff = 200; end
    if k < beats.count
        tNext = beats.times(k+1);
    else
        tNext = res.t(end);
    end

    try
        idxAt = @(tt) find(res.t >= min(tt, res.t(end)), 1, 'first');
        pts = [idxAt(tCross), idxAt(tPeak + 0.3*apdEff), idxAt(tCross + apdEff), ...
               idxAt((tCross + apdEff + tNext)/2)];
        labels = {'0','2','3','4'};
        for j = 1:numel(pts)
            plot(res.t(pts(j)), res.V(pts(j)), 'o', 'MarkerFaceColor',[0.85 0.1 0.1], ...
                'MarkerEdgeColor','k', 'MarkerSize',6);
            text(res.t(pts(j)), res.V(pts(j))+6, labels{j}, 'FontWeight','bold', ...
                'FontSize',12, 'HorizontalAlignment','center');
        end
        annotation('textbox', [0.58 0.14 0.40 0.32], 'String', { ...
            '0: rapid upstroke (fast Na^+ current)', ...
            '2: plateau (persistent Na^+ + slowly rising i_{K2})', ...
            '3: repolarization (i_{K2}, i_{K1} dominate)', ...
            '4: diastolic depolarization (decaying i_{K2})', ...
            '(no phase-1 notch: this 1962 model predates I_{to})'}, ...
            'EdgeColor','none', 'FontSize',8.5, 'BackgroundColor','w');
    catch
    end
end
end
