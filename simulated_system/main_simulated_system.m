clear
clc

h = 0.01;
Tsim = 10;
x0 = [0 0 0 0]';

t_u = (0:h:Tsim)';
u = timeseries(1*sin(2*pi*t_u), t_u);

param.M   = 0.2;     % cart mass [kg]
param.m   = 0.05;    % pendulum mass [kg]
param.g   = 9.81;    % gravity [m/s^2]
param.b   = 0.001;   % pendulum damping
param.c   = 0.005;   % cart damping
param.l   = 0.2;     % pendulum length [m]
param.k_m = 1.0;     % motor gain

disp('simulating ...')

output = sim('simulated_system');

disp('plotting ...')

visualize_trajectory(output.tout, output.x, param)

disp('saving data ...')

results_dir = fullfile(pwd, 'sim_results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir)
end

results_dir = fullfile(pwd, 'sim_results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir)
end

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
filename = sprintf('data_%s.mat', timestamp);
save(fullfile(results_dir, filename), 'output', 'u', 'param')

disp('DONE')