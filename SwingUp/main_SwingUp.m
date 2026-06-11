h = 0.01;
T = 10;


timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

% add path to the EOM and other functions, assuming this script is in LQI/
proj_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(proj_root, 'global'));

results_dir = fullfile(proj_root, 'SwingUp', 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
run_dir = fullfile(results_dir, timestamp);
if ~exist(run_dir, 'dir'), mkdir(run_dir); end


load("global/swing_up.mat");
load("global/param_64_175.mat");
fields = fieldnames(param);
for i = 1:numel(fields)
    assignin("base",fields{i}, param.(fields{i}))
end

mname = 'inverted_pendulum_SU';
play_run(mname);

fprintf('Waiting for experiment to finish...\n');
while ~strcmp(get_param(mname, 'SimulationStatus'), 'stopped')
    pause(1);
end
fprintf('Experiment complete.\n');



% Collect workspace variables written by the To Workspace blocks
results_hw.t          = 0:h:T;
results_hw.x          = evalin('base', 'x');
results_hw.theta      = evalin('base', 'theta');
results_hw.x_hat      = evalin('base', 'x_hat');
results_hw.theta_hat  = evalin('base', 'theta_hat');
results_hw.dx_hat     = evalin('base', 'dx_hat');
results_hw.dtheta_hat = evalin('base', 'dtheta_hat');

res_file = fullfile(run_dir, ['run_hw_' timestamp '.mat']);
save(res_file, '-struct', 'results_hw');
fprintf('Hardware results saved to %s\n', res_file);


plot_kalman_comparison(results_hw.x, results_hw.theta, ...
    results_hw.x_hat, results_hw.theta_hat, ...
    results_hw.dx_hat, results_hw.dtheta_hat, ...
    h, T, run_dir, timestamp);
plot_augmented_states(results_hw.x_hat, results_hw.dx_hat, ...
    results_hw.theta_hat, results_hw.dtheta_hat, ...
    h, T, run_dir, timestamp);