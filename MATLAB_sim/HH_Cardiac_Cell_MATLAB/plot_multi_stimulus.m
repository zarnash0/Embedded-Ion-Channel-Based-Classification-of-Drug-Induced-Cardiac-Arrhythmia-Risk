function fig = plot_multi_stimulus(params, Iamps, opts)
arguments
    params (1,1) struct
    Iamps (1,:) double
    opts.tEnd (1,1) double = 100
    opts.onset (1,1) double = 10
    opts.offset (1,1) double = 90
    opts.dt (1,1) double = 0.01
end

fig = figure('Color','w','Position',[100 100 950 560]);
cmap = turbo(numel(Iamps));
hold on;
for k = 1:numel(Iamps)
    stim = make_pulse(opts.onset, opts.offset, Iamps(k));
    res  = simulate_HH(params, stim, opts.tEnd, dt=opts.dt);
    plot(res.t, res.V, 'LineWidth', 1.3, 'Color', cmap(k,:), ...
         'DisplayName', sprintf('%.0f nA/mm^2', Iamps(k)));
end
xlabel('Time (ms)'); ylabel('Membrane potential V_m (mV)');
title('Hodgkin-Huxley response across stimulus amplitudes');
legend('Location','eastoutside'); grid on; box off;
end
