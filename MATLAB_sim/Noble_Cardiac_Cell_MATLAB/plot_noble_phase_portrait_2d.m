function fig = plot_noble_phase_portrait_2d(res)
pal = noble_palette();
dVdt = [0, diff(res.V)./diff(res.t)];
depolMask = dVdt >= 0;

fig = figure('Color','w','Position',[100 100 750 600]);
scatter(res.V(depolMask), res.n(depolMask), 4, 'filled', ...
    'MarkerFaceColor', pal.depol, 'MarkerFaceAlpha', 0.6); hold on;
scatter(res.V(~depolMask), res.n(~depolMask), 4, 'filled', ...
    'MarkerFaceColor', pal.repol, 'MarkerFaceAlpha', 0.6);
xlabel('V_m (mV)'); ylabel('n (K^+ activation, i_{K2})');
title('Phase portrait: V - n, colored by dV/dt sign (Fenton/Virtual-Heart convention)');
legend({'depolarizing (dV/dt \geq 0)','repolarizing (dV/dt < 0)'}, 'Location','best');
grid on; box off;
end
