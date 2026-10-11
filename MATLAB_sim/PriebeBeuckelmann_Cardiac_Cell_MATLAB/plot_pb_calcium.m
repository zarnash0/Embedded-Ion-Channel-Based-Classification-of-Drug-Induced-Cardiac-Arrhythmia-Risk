function fig = plot_pb_calcium(res)
%       1. V_m, for temporal reference.
%       2. [Ca^{2+}]_i in nM, with its diastolic and systolic values
%          annotated on the trace.
%       3. the two SR compartments, [Ca^{2+}]_{JSR} and [Ca^{2+}]_{NSR}.
%       4. the SR fluxes I_up (SERCA uptake), I_leak and I_tr, together
%          with the two sarcolemmal currents that actually carry Ca2+
%          across the membrane, I_Ca and I_NaCa (two y-axes: the fluxes
%          and the currents are in genuinely different units).
arguments
    res (1,1) struct
end

if ~isfield(res, 'model') || ~strcmp(string(res.model), "PB_full")
    if isfield(res, 'model'), got = char(string(res.model)); else, got = '<no .model field>'; end
    error('plot_pb_calcium:wrongModel', ...
        ['PLOT_PB_CALCIUM applies only to the full model (res.model == ''PB_full''); got ''%s''.\n' ...
         'The Bernus et al. (2002) reduced model has no calcium handling to plot: it FIXES\n' ...
         '[Ca2+]i at a constant 400 nM (and [Na+]i / [K+]i at 10 / 140 mM) and removes the\n' ...
         'JSR/NSR compartments and their fluxes entirely. That removal is the defining\n' ...
         'simplification which takes the model from 15 state variables to 6, so there is no\n' ...
         '[Ca2+]i transient, no CaJSR/CaNSR and no Iup/Ileak/Itr in its output struct.\n' ...
         'Use PLOT_PB_PHASE3D instead: it branches on res.model and plots (V, X, f) for the\n' ...
         'reduced model, i.e. the two slow gates that shape its plateau.'], got);
end

pal = pb_palette();
t      = res.t;
Cai_nM = res.Cai * 1e6;          % model stores [Ca2+]i in mM; 1 mM = 1e6 nM

beats = detect_pb_ap(res);
[iA, iB, kBeat, CL] = local_beat_window(t, beats);
w  = iA:iB;
tw = t(w);

if isfield(res, 'params') && isfield(res.params, 'condition')
    condLabel = char(string(res.params.condition));
else
    condLabel = 'unspecified condition';
end

fig = figure('Color','w','Position',[80 40 1000 940]);
tl = tiledlayout(5, 1, 'TileSpacing','compact');
title(tl, sprintf('Ca^{2+} handling - full Priebe-Beuckelmann (1998) model (%s)', condLabel), ...
    'FontWeight','bold');

% ---------------- Panel 1: membrane potential --------------------------
ax1 = nexttile;
plot(tw, res.V(w), 'Color', pal.V, 'LineWidth', 1.6);
hold on;
if ~isnan(kBeat)
    % Stimulus onset, so the reader can line up every panel below against
    % the same t = 0 of the beat.
    xline(beats.stimTimes(kBeat), '-', 'Color', pal.stim, 'LineWidth', 1.0);
end
ylabel('V_m (mV)');
title('Membrane potential (last complete beat)');
grid on; box off; xlim([tw(1) tw(end)]);

% ---------------- Panel 2: cytosolic calcium transient -----------------
ax2 = nexttile;
plot(tw, Cai_nM(w), 'Color', pal.ICa, 'LineWidth', 1.6);
hold on;

% Diastolic value is taken at the stimulus instant (the true pre-upstroke
% baseline of THIS beat), not as min() over the window, which would land
% late in the previous beat's decay if the cell has not fully relaxed.
if ~isnan(kBeat)
    iStim = find(t >= beats.stimTimes(kBeat), 1, 'first');
else
    iStim = iA;
end
CaiRest = Cai_nM(iStim);
[CaiPeak, iPkRel] = max(Cai_nM(w));
tPk = tw(iPkRel);

% REFERENCE VALUES (deliberately NOT drawn on the figure, so the plot
% always shows what this implementation actually produced): Priebe &
% Beuckelmann report ~120 nM diastolic / 614 nM peak for the nonfailing
% cell and 136 nM / 334 nM for the failing cell at 1 Hz steady state.
% Our implementation gives ~709 nM peak nonfailing - same shape and
% amplitude class, so the comparison belongs in the write-up, not baked
% into an annotation that could go stale.
yline(CaiRest, ':', 'Color', pal.grid, 'LineWidth', 1.0);
plot(tPk, CaiPeak, 'o', 'MarkerFaceColor', pal.ICa, 'MarkerEdgeColor', pal.V, 'MarkerSize', 6);
plot(t(iStim), CaiRest, 's', 'MarkerFaceColor', pal.grid, 'MarkerEdgeColor', pal.V, 'MarkerSize', 6);
text(tPk, CaiPeak, sprintf('  peak %.0f nM', CaiPeak), ...
    'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize', 9);
text(tw(1), CaiRest, sprintf('  diastolic %.0f nM  ', CaiRest), ...
    'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize', 9);
text(0.99, 0.92, sprintf('amplitude %.0f nM', CaiPeak - CaiRest), 'Units','normalized', ...
    'HorizontalAlignment','right', 'FontSize', 9, 'BackgroundColor','w');

ylabel('[Ca^{2+}]_i (nM)');
title('Cytosolic Ca^{2+} transient');
grid on; box off; xlim([tw(1) tw(end)]);

% ---------------- Panel 3: SR calcium stores ---------------------------
ax3 = nexttile;
hJ = plot(tw, res.CaJSR(w), 'Color', pal.to, 'LineWidth', 1.5, 'LineStyle','-'); hold on;
hN = plot(tw, res.CaNSR(w), 'Color', pal.X,  'LineWidth', 1.5, 'LineStyle','-');
ylabel('SR [Ca^{2+}] (mM)');
% The JSR dip immediately after the upstroke IS the CICR release event:
% SIMULATE_PB_FULL does not export I_rel as a field, so the release flux
% is read here indirectly, as junctional depletion followed by refilling
% from the network SR through I_tr.
title('SR stores: junctional (JSR) and network (NSR) - the JSR dip is the CICR release');
legend([hJ hN], {'[Ca^{2+}]_{JSR}','[Ca^{2+}]_{NSR}'}, 'Location','best', 'Box','off');
grid on; box off; xlim([tw(1) tw(end)]);

% ---------------- Panel 4: calcium fluxes and currents -----------------
% Two y-axes because these are not the same physical quantity: I_up,
% I_leak and I_tr are compartmental Ca2+ fluxes in mM/ms (shown as uM/ms
% for readability, i.e. x1e3), whereas I_Ca and I_NaCa are transmembrane
% current densities in pA/pF. Plotting them on a shared axis would be
% meaningless. Note also that these are the raw flux expressions from
% PB_FULL_CURRENTS, before the volume ratios (V_JSR/V_myo, V_NSR/V_myo)
% that scale them in the concentration ODEs.
ax4 = nexttile;
yyaxis left
hUp   = plot(tw, res.Iup(w)*1e3,   'Color', pal.X,    'LineWidth', 1.4, 'LineStyle','-'); hold on;
hLeak = plot(tw, res.Ileak(w)*1e3, 'Color', pal.m,    'LineWidth', 1.4, 'LineStyle','-');
hTr   = plot(tw, res.Itr(w)*1e3,   'Color', pal.to,   'LineWidth', 1.4, 'LineStyle','-');
ylabel('SR flux (\muM/ms)');
yyaxis right
hICa  = plot(tw, res.ICa(w),       'Color', pal.ICa,  'LineWidth', 1.4, 'LineStyle','-');
hNaCa = plot(tw, res.INaCa(w),     'Color', pal.INaCa,'LineWidth', 1.4, 'LineStyle','-');
ylabel('Current (pA/pF)');
% Fluxes are coloured by the compartment they feed (blue = NSR, matching
% CaNSR above; dark red = JSR, matching CaJSR), and the currents keep the
% project-wide pb_palette identities for I_Ca and I_NaCa.
ax4.YAxis(1).Color = pal.V;
ax4.YAxis(2).Color = pal.V;
legend([hUp hLeak hTr hICa hNaCa], ...
    {'I_{up} (SERCA)','I_{leak}','I_{tr} (NSR\rightarrowJSR)','I_{Ca}','I_{NaCa}'}, ...
    'Location','best', 'Box','off', 'NumColumns', 2);
xlabel('Time (ms)');
title('Ca^{2+} fluxes (left axis) and Ca^{2+}-carrying membrane currents (right axis)');
grid on; box off; xlim([tw(1) tw(end)]);

% ---------------- Panel 5: whole-run drift -----------------------------
ax5 = nexttile;
% The trace is decimated for rendering speed only; the beat-by-beat
% envelope below is computed from the FULL-resolution data, so no peak
% can be missed by the decimation.
stride = max(1, ceil(numel(t)/20000));
hTrace = plot(t(1:stride:end), Cai_nM(1:stride:end), 'Color', pal.ICa, 'LineWidth', 1.0);
hold on;
[dias, sys] = local_beat_envelope(t, Cai_nM, beats, CL);
handles = hTrace;
labels  = {'[Ca^{2+}]_i'};
if any(isfinite(sys))
    hSys = plot(beats.stimTimes, sys, 'o-', 'Color', pal.full, 'MarkerFaceColor', pal.full, ...
        'MarkerSize', 4, 'LineWidth', 1.2);
    handles(end+1) = hSys; labels{end+1} = 'per-beat peak';
end
if any(isfinite(dias))
    hDia = plot(beats.stimTimes, dias, 's-', 'Color', pal.reduced, 'MarkerFaceColor', pal.reduced, ...
        'MarkerSize', 4, 'LineWidth', 1.2);
    handles(end+1) = hDia; labels{end+1} = 'per-beat diastolic';
end

yl = ylim(ax5);
% Shade the beat that panels 1-4 zoom in on, so the zoom is locatable
% within the whole record.
hPatch = patch([t(iA) t(iB) t(iB) t(iA)], [yl(1) yl(1) yl(2) yl(2)], pal.grid, ...
    'EdgeColor','none', 'FaceAlpha', 0.6);
uistack(hPatch, 'bottom');
ylim(ax5, yl);

if beats.count >= 2
    kA = find(isfinite(sys), 1, 'first');
    kB = find(isfinite(sys), 1, 'last');
    if ~isempty(kA) && ~isempty(kB) && kB > kA
        text(0.015, 0.90, sprintf(['drift over %d beats:  peak %.0f \\rightarrow %.0f nM,' ...
            '  diastolic %.0f \\rightarrow %.0f nM'], kB-kA+1, sys(kA), sys(kB), dias(kA), dias(kB)), ...
            'Units','normalized', 'FontSize', 9, 'BackgroundColor','w');
    end
end

legend(handles, labels, 'Location','best', 'Box','off', 'NumColumns', 3);
xlabel('Time (ms)'); ylabel('[Ca^{2+}]_i (nM)');
title('Whole run: beat-to-beat Ca^{2+} drift (a key motivation for the reduced model)');
grid on; box off; xlim([t(1) t(end)]);

linkaxes([ax1 ax2 ax3 ax4], 'x');
end

% =======================================================================
function [iA, iB, kBeat, CL] = local_beat_window(t, beats)
%LOCAL_BEAT_WINDOW Index range covering the last COMPLETE paced beat.
%   Falls back to the whole record if no stimuli were detected.
if beats.count == 0
    iA = 1; iB = numel(t); kBeat = NaN; CL = t(end) - t(1);
    return
end
if beats.count >= 2
    CL = median(diff(beats.stimTimes));
else
    CL = max(t(end) - beats.stimTimes(1), eps);
end
% Prefer the final beat, but step back one if the record ends before that
% beat has run a full cycle: otherwise the panels would show a
% repolarization and a Ca2+ decay chopped off by the end of the
% simulation, which reads as a modelling artefact rather than a plotting
% one.
kBeat = beats.count;
if beats.count >= 2 && (t(end) - beats.stimTimes(end)) < 0.95*CL
    kBeat = beats.count - 1;
end
tA = max(t(1),   beats.stimTimes(kBeat) - 0.05*CL);   % a little pre-stimulus baseline
tB = min(t(end), beats.stimTimes(kBeat) + CL);
iA = find(t >= tA, 1, 'first');
iB = find(t <= tB, 1, 'last');
if isempty(iA), iA = 1; end
if isempty(iB) || iB <= iA, iB = numel(t); end
end

% =======================================================================
function [dias, sys] = local_beat_envelope(t, Cai_nM, beats, CL)
%LOCAL_BEAT_ENVELOPE Diastolic and systolic [Ca2+]i for every paced beat.
%   Computed on the full-resolution trace (see the decimation note above).
dias = nan(1, beats.count);
sys  = nan(1, beats.count);
for k = 1:beats.count
    i0 = find(t >= beats.stimTimes(k), 1, 'first');
    if k < beats.count
        i1 = find(t < beats.stimTimes(k+1), 1, 'last');
    else
        i1 = find(t <= beats.stimTimes(k) + CL, 1, 'last');
    end
    if isempty(i0) || isempty(i1) || i1 <= i0
        continue
    end
    dias(k) = Cai_nM(i0);          % value at the stimulus = pre-upstroke baseline
    sys(k)  = max(Cai_nM(i0:i1));
end
end
