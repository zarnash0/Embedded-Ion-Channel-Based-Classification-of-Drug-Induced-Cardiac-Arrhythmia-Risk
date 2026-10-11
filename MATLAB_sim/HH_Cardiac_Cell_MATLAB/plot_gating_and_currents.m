function fig = plot_gating_and_currents(res)
pal = hh_palette();
fig = figure('Color','w','Position',[100 100 900 620]);
tiledlayout(2,1,'TileSpacing','compact');

nexttile
plot(res.t, res.m, 'Color', pal.m, 'LineWidth', 1.5); hold on;
plot(res.t, res.h, 'Color', pal.h, 'LineWidth', 1.5);
plot(res.t, res.n, 'Color', pal.n, 'LineWidth', 1.5);
ylim([0 1]); ylabel('Gating variable');
title('Gating variables m, h, n');
legend({'m (Na^+ activation)','h (Na^+ inactivation)','n (K^+ activation)'}, 'Location','best');
grid on; box off;

nexttile
plot(res.t, res.INa, 'Color', pal.INa, 'LineWidth', 1.3); hold on;
plot(res.t, res.IK,  'Color', pal.IK,  'LineWidth', 1.3);
plot(res.t, res.IL,  'Color', pal.IL,  'LineWidth', 1.3);
plot(res.t, res.INa+res.IK+res.IL, 'k--', 'LineWidth', 1);
xlabel('Time (ms)'); ylabel('Ionic current (nA/mm^2)');
title('Ionic membrane currents');
legend({'I_{Na}','I_K','I_L','I_{ion} total'}, 'Location','best');
grid on; box off;
end
