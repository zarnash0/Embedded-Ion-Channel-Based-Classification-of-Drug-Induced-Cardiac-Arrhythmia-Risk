function fig = plot_phase_portraits(res)
pal = hh_palette();
fig = figure('Color','w','Position',[100 100 950 460]);
tiledlayout(1,2,'TileSpacing','compact');

nexttile
plot(res.V, res.n, 'Color', pal.n, 'LineWidth', 1.4);
xlabel('V_m (mV)'); ylabel('n');
title('Phase portrait: V - n (K^+ activation)');
grid on; box off;

nexttile
plot(res.V, res.m.^3 .* res.h, 'Color', pal.m, 'LineWidth', 1.4);
xlabel('V_m (mV)'); ylabel('m^3h');
title('Phase portrait: V - m^3h (Na^+ availability)');
grid on; box off;
end
