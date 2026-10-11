function fig = plot_FI_curve(fi)
%PLOT_FI_CURVE Sustained firing rate vs. injected current.
pal = hh_palette();
fig = figure('Color','w','Position',[100 100 700 500]);
plot(fi.I, fi.rate, 'o-', 'Color', pal.n, 'LineWidth', 1.6, ...
     'MarkerFaceColor', pal.n, 'MarkerSize', 5);
xlabel('Injected current I_e/A (nA/mm^2)');
ylabel('Sustained firing rate (Hz)');
title('F-I curve');
grid on; box off;
end
