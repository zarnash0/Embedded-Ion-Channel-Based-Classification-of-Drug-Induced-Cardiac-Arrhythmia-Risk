function fig = plot_noble_gating_and_currents(res)
pal = noble_palette();
fig = figure('Color','w','Position',[100 100 950 680]);
tiledlayout(2,1,'TileSpacing','compact');

nexttile
plot(res.t, res.m, 'Color', pal.m, 'LineWidth', 1.4); hold on;
plot(res.t, res.h, 'Color', pal.h, 'LineWidth', 1.4);
plot(res.t, res.n, 'Color', pal.n, 'LineWidth', 1.4);
ylim([0 1]); ylabel('Gating variable');
title('Gating variables m, h, n');
legend({'m (Na^+ activation)','h (Na^+ inactivation)','n (K^+ activation, i_{K2})'}, 'Location','best');
grid on; box off;

nexttile
plot(res.t, res.iNa, 'Color', pal.iNa, 'LineWidth', 1.2); hold on;
plot(res.t, res.iK1, 'Color', pal.iK1, 'LineWidth', 1.2);
plot(res.t, res.iK2, 'Color', pal.iK2, 'LineWidth', 1.2);
plot(res.t, res.iAn, 'Color', pal.iAn, 'LineWidth', 1.2);
plot(res.t, res.iNa+res.iK1+res.iK2+res.iAn, 'k--', 'LineWidth', 1);
xlabel('Time (ms)'); ylabel('Ionic current (nA)');
title('Ionic membrane currents');
legend({'i_{Na}','i_{K1} (anomalous rectifier)','i_{K2} (slow, gated by n)','i_{An} (leak)','i_{ion} total'}, 'Location','best');
grid on; box off;
end
