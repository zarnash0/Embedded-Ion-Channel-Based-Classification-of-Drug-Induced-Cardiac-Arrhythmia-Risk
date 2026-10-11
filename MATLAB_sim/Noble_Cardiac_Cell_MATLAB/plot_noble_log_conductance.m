function fig = plot_noble_log_conductance(res, opts)
%   Replicates Fink & Noble (2008, Scholarpedia "Noble model") Figure 4:
%   "Sodium and potassium conductance changes computed from the model,"
%   plotted on a log axis to accommodate the large changes in sodium
%   conductance. This is one of the two figures the source article
%   itself uses as validation
arguments
    res (1,1) struct
    opts.nCycles (1,1) double = 2
    opts.threshold (1,1) double = -40
end

pal = noble_palette();
beats = detect_beats(res, threshold=opts.threshold);
if beats.count >= opts.nCycles+1
    tStart = beats.times(1);
    tEnd   = beats.times(min(beats.count, 1+opts.nCycles));
else
    tStart = res.t(1);
    tEnd   = res.t(end);
end
idxRange = res.t >= tStart & res.t <= tEnd;

gK = res.gK1(idxRange) + res.gK2(idxRange);
allVals = [res.gNa(idxRange), gK];

fig = figure('Color','w','Position',[100 100 950 500]);
semilogy(res.t(idxRange), res.gNa(idxRange), 'Color', pal.iNa, 'LineWidth', 1.4); hold on;
semilogy(res.t(idxRange), gK, 'Color', pal.iK2, 'LineWidth', 1.4);
xlabel('Time (ms)'); ylabel('Conductance (\muS), log scale');
title('Sodium and potassium conductance changes (cf. Fink & Noble 2008, Fig. 4)');
legend({'g_{Na} (= 400000\cdotm^3h + 140)','g_K (= g_{K1}+g_{K2})'}, 'Location','best');

ylim([0.7*min(allVals(allVals>0)), 1.3*max(allVals)]);
grid on; box off;
end
