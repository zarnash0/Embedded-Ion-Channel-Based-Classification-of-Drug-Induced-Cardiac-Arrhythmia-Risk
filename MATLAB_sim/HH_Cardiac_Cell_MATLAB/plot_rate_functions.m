function fig = plot_rate_functions()
pal = hh_palette();
Vr = linspace(-100, 50, 400);
[a, b] = hh_gating_rates(Vr);
minf = a.m./(a.m+b.m); taum = 1./(a.m+b.m);
hinf = a.h./(a.h+b.h); tauh = 1./(a.h+b.h);
ninf = a.n./(a.n+b.n); taun = 1./(a.n+b.n);

fig = figure('Color','w','Position',[100 100 950 680]);
tiledlayout(2,2,'TileSpacing','compact');

nexttile
plot(Vr, minf, 'Color', pal.m, 'LineWidth', 1.6); hold on;
plot(Vr, hinf, 'Color', pal.h, 'LineWidth', 1.6);
plot(Vr, ninf, 'Color', pal.n, 'LineWidth', 1.6);
ylabel('x_\infty'); title('Steady-state (in)activation');
legend({'m_\infty','h_\infty','n_\infty'}, 'Location','best'); grid on; box off;

nexttile
plot(Vr, taum, 'Color', pal.m, 'LineWidth', 1.6); hold on;
plot(Vr, tauh, 'Color', pal.h, 'LineWidth', 1.6);
plot(Vr, taun, 'Color', pal.n, 'LineWidth', 1.6);
ylabel('\tau (ms)'); title('Time constants');
legend({'\tau_m','\tau_h','\tau_n'}, 'Location','best'); grid on; box off;

nexttile
plot(Vr, a.m, 'Color', pal.m, 'LineWidth', 1.4); hold on;
plot(Vr, a.h, 'Color', pal.h, 'LineWidth', 1.4);
plot(Vr, a.n, 'Color', pal.n, 'LineWidth', 1.4);
xlabel('V_m (mV)'); ylabel('\alpha (1/ms)'); title('\alpha rate constants');
legend({'\alpha_m','\alpha_h','\alpha_n'}, 'Location','best'); grid on; box off;

nexttile
plot(Vr, b.m, 'Color', pal.m, 'LineWidth', 1.4); hold on;
plot(Vr, b.h, 'Color', pal.h, 'LineWidth', 1.4);
plot(Vr, b.n, 'Color', pal.n, 'LineWidth', 1.4);
xlabel('V_m (mV)'); ylabel('\beta (1/ms)'); title('\beta rate constants');
legend({'\beta_m','\beta_h','\beta_n'}, 'Location','best'); grid on; box off;
end
