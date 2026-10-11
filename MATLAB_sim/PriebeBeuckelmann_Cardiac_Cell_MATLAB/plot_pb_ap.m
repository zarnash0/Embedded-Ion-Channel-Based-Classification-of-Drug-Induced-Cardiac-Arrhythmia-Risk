function fig = plot_pb_ap(res, opts)
%     beatIndex      which beat to annotate. 0 (default) selects the last
%                    beat whose repolarization is fully contained in the
%                    record -- normally the closest thing to steady state
%                    that the simulation offers.
%     repolFraction  repolarization fraction used for the bracket
%                    [default 0.90, i.e. the usual APD90]
arguments
    res (1,1) struct
    opts.beatIndex (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    opts.repolFraction (1,1) double {mustBePositive, mustBeLessThan(opts.repolFraction,1)} = 0.90
end

pal = pb_palette();
t = res.t;
V = res.V;
ap = detect_pb_ap(res, repolFractions=opts.repolFraction);

k = 0;
if ap.count > 0
    if opts.beatIndex == 0
        good = find(isfinite(ap.APD(:,1)));
        if isempty(good), k = ap.count; else, k = good(end); end
    elseif opts.beatIndex > ap.count
        warning('plot_pb_ap:beatIndex', ...
            'beatIndex %d exceeds the %d detected beats; annotating the last one instead.', ...
            opts.beatIndex, ap.count);
        k = ap.count;
    else
        k = opts.beatIndex;
    end
end

Vmin = min(V); Vmax = max(V);
span = max(Vmax - Vmin, 1);          % guard against a perfectly flat trace
bandBase = Vmin - 0.10*span;         % zero level of the display stimulus trace
bandH    = 0.055*span;               % height of a full-amplitude stimulus pulse
yLo      = bandBase - 0.05*span;
yHi      = Vmax + 0.11*span;

fig = figure('Color','w','Position',[100 100 980 560]);
hold on;
if ap.count > 0 && ap.count <= 30
    for i = 1:ap.count
        xline(ap.stimTimes(i), ':', 'Color', pal.stim, 'LineWidth', 0.5, ...
              'HandleVisibility','off');
    end
end
Iref = max(abs(res.Istim));
if Iref > 0
    stimTrace = bandBase + bandH*(res.Istim/Iref);
    plot(t, stimTrace, '-', 'Color', pal.stim, 'LineWidth', 1.1, ...
        'DisplayName', sprintf('I_{stim} (peak %.0f pA/pF, offset for display)', Iref));
end

% The action potential
plot(t, V, '-', 'Color', pal.V, 'LineWidth', 1.6, 'DisplayName', 'V_m');

% Beat annotation
if k > 0
    RMP  = ap.RMP(k);
    OS   = ap.overshoot(k);
    amp  = ap.amplitude(k);
    tStim = ap.stimTimes(k);
    yline(OS, ':', sprintf('overshoot %+.1f mV', OS), 'Color', pal.stim, ...
        'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top', ...
        'FontSize', 8.5, 'HandleVisibility','off');
    yline(RMP, ':', sprintf('RMP %.1f mV', RMP), 'Color', pal.stim, ...
        'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','bottom', ...
        'FontSize', 8.5, 'HandleVisibility','off');

    if isfinite(ap.APD(k,1))
        apd  = ap.APD(k,1);
        yLvl = OS - opts.repolFraction*amp;      % the repolarization level 
        tick = 0.022*span;

        plot([tStim tStim+apd], [yLvl yLvl], '-', 'Color', pal.V, ...
             'LineWidth', 1.1, 'HandleVisibility','off');
        plot([tStim tStim], yLvl+[-tick tick], '-', 'Color', pal.V, ...
             'LineWidth', 1.1, 'HandleVisibility','off');
        plot([tStim+apd tStim+apd], yLvl+[-tick tick], '-', 'Color', pal.V, ...
             'LineWidth', 1.1, 'HandleVisibility','off');
        text(tStim+apd/2, yLvl+2.4*tick, ...
            sprintf('APD_{%d} = %.1f ms', round(100*opts.repolFraction), apd), ...
            'HorizontalAlignment','center', 'FontSize', 9, 'BackgroundColor','w');
    end
end

xlim([t(1) t(end)]);
ylim([yLo yHi]);
xlabel('Time (ms)');
ylabel('Membrane potential V_m (mV)');
grid on; box off;
legend('Location','best', 'Box','off');

if k > 0
    if isfinite(ap.APD(k,1))
        apdStr = sprintf('APD_{%d} %.1f ms', round(100*opts.repolFraction), ap.APD(k,1));
    else
        apdStr = 'APD not measurable (record ends before repolarization)';
    end
    sub = sprintf('beat %d of %d:  RMP %.1f mV,  overshoot %+.1f mV,  amplitude %.1f mV,  dV/dt_{max} %.0f V/s,  %s', ...
        k, ap.count, ap.RMP(k), ap.overshoot(k), ap.amplitude(k), ap.dVdtMax(k), apdStr);
else
    sub = 'no paced beats detected in this record';
end
title({localModelLabel(res), sub});
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
