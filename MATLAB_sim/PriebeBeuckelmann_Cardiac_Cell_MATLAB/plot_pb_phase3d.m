function fig = plot_pb_phase3d(res, opts)
%PLOT_PB_PHASE3D 3D single-cell state-space trajectory, coloured by time.
%   FIG = PLOT_PB_PHASE3D(RES) plots the trajectory of the simulation
%   result RES through a three-dimensional slice of its own state space,
%   with the scatter coloured by simulation time so the direction of
%   travel and the beat-to-beat repetition are both readable. The axes
%   are chosen per model (branching on RES.model):
%
%     'PB_full'    (V, [Ca^{2+}]_i, [Ca^{2+}]_{JSR}). These three carry
%                  the calcium-induced calcium release loop directly:
%                  depolarization opens I_Ca, the resulting Ca2+ influx
%                  triggers release from the junctional SR (JSR depletes),
%                  cytosolic Ca2+ rises, and SERCA refills the store
%                  during diastole. The loop closes visibly in 3D in a
%                  way it does not in any pair of 2D projections.
%     'PB_reduced' (V, X, f). The reduced model fixes [Ca^{2+}]_i at a
%                  constant 400 nM, so calcium axes are degenerate;
%                  instead this plots the two SLOW gates that shape its
%                  plateau - X, the combined delayed-rectifier activation
%                  gate, and f, the L-type Ca2+ inactivation gate.
%
%   Name-value options:
%     maxPoints     maximum number of plotted samples        [default 20000]
%     colormapName  "turbo" (default) or "parula"
%
%   ORIGINALITY NOTE (stated plainly, as elsewhere in this project): a 3D
%   single-cell state-space trajectory of this kind is an ORIGINAL
%   presentation choice for this project. It is NOT a reproduction of any
%   figure in Priebe & Beuckelmann (1998) or in Bernus et al. (2002) -
%   neither paper contains a phase-space rendering of any kind - and it
%   is labelled as such on the figure itself. It follows the same
%   convention already established at the Noble stage of this project
%   (PLOT_NOBLE_PHASE_PORTRAIT_3D), so that all three models' figures
%   read as one visual system.
%
%   Priebe L & Beuckelmann DJ (1998). Circ Res 82:1206-1223.
%   Bernus O, Wilders R, Zemlin CW, Verschelde H, Panfilov AV (2002).
%   Am J Physiol Heart Circ Physiol 282:H2296-H2308.
%
%   See also SIMULATE_PB_FULL, SIMULATE_PB_REDUCED, PLOT_PB_CALCIUM.

arguments
    res (1,1) struct
    opts.maxPoints (1,1) double {mustBePositive} = 20000
    opts.colormapName (1,1) string {mustBeMember(opts.colormapName, ["turbo","parula"])} = "turbo"
end

if ~isfield(res, 'model')
    error('plot_pb_phase3d:noModel', ...
        'RES has no .model field; expected the output of SIMULATE_PB_FULL or SIMULATE_PB_REDUCED.');
end
model = string(res.model);

switch model
    case "PB_full"
        xv = res.V;
        yv = res.Cai * 1e6;          % mM -> nM, for a readable axis scale
        zv = res.CaJSR;              % mM
        xlab = 'V_m (mV)';
        ylab = '[Ca^{2+}]_i (nM)';
        zlab = '[Ca^{2+}]_{JSR} (mM)';
        if isfield(res, 'params') && isfield(res.params, 'condition')
            variant = char(string(res.params.condition));
        else
            variant = 'unspecified condition';
        end
        headline = sprintf('Full Priebe-Beuckelmann (1998), %s: state space (V, [Ca^{2+}]_i, [Ca^{2+}]_{JSR})', variant);
        subline  = 'the calcium-induced calcium release loop, traced directly';

    case "PB_reduced"
        % [Ca2+]i is a constant in this model, so a calcium axis would
        % collapse to a plane; the two slow gates are the informative
        % pair here.
        xv = res.V;
        yv = res.X;
        zv = res.f;
        xlab = 'V_m (mV)';
        ylab = 'X  (delayed-rectifier activation)';
        zlab = 'f  (I_{Ca} inactivation)';
        if isfield(res, 'params') && isfield(res.params, 'cellType')
            variant = char(string(res.params.cellType));
        else
            variant = 'unspecified cell type';
        end
        headline = sprintf('Reduced Bernus et al. (2002), %s: state space (V, X, f)', variant);
        subline  = '[Ca^{2+}]_i is fixed at 400 nM in this model, so the slow gates replace it';

    otherwise
        error('plot_pb_phase3d:unknownModel', ...
            ['RES.model is ''%s''; expected ''PB_full'' or ''PB_reduced''. ' ...
             'This function only understands the two Priebe-Beuckelmann-family models ' ...
             'implemented in this project.'], char(model));
end

% Uniform decimation for rendering speed. This is a display-only measure:
% a phase-space trajectory is a geometric object, so dropping intermediate
% samples thins the point cloud without moving the curve, unlike a time
% series where a decimated peak reads as a smaller peak.
n = numel(res.t);
stride = max(1, ceil(n/opts.maxPoints));
idx = 1:stride:n;

fig = figure('Color','w','Position',[100 100 900 700]);
scatter3(xv(idx), yv(idx), zv(idx), 5, res.t(idx), 'filled');
colormap(fig, char(opts.colormapName));
cb = colorbar; cb.Label.String = 'Time (ms)';
xlabel(xlab); ylabel(ylab); zlabel(zlab);
title({headline, subline, ...
       '(original presentation choice for this project - not a figure from either source paper)'});
grid on; box off;
view(-35, 20);
end
