function [out] = simulate_system(h, Tsim, x0, u, param, save_results, view_traj)
simIn = Simulink.SimulationInput('simulated_system');

simIn = simIn.setModelParameter('StopTime',  num2str(Tsim));
simIn = simIn.setModelParameter('FixedStep', num2str(h));

simIn = simIn.setVariable('x0',     x0);
simIn = simIn.setVariable('u',      u);
simIn = simIn.setVariable('params', param);

out = sim(simIn);

if nargin <= 6
    view_traj = false;
    if naargin <= 5
        save_results = false;
    end
end

if save_results
    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    filename = sprintf('sim_%s.mat', timestamp);
    save(fullfile(pwd, 'simulated_system', 'sim_results', filename), ...
        'out', 'u', 'param')
end

if view_traj
    visualize_trajectory(out.tout, out.x, param)
end
