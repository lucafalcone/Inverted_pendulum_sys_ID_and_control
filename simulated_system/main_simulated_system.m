clear
clc

h = 0.01;
Tsim = 10;
x0 = [0 0 0 0]';

t_u = (0:h:Tsim)';
u = timeseries(1*sin(1*pi*t_u), t_u);

param.M   = 0.2;     % cart mass [kg]
param.m   = 0.05;    % pendulum mass [kg]
param.g   = 9.81;    % gravity [m/s^2]
param.b   = 0.001;   % pendulum damping
param.c   = 0.005;   % cart damping
param.l   = 0.2;     % pendulum length [m]
param.k_m = 1.0;     % motor gain

save_results = true;
view_traj = true;
out = simulate_system(h, Tsim, x0, u, param, save_results, view_traj);


% disp('simulating ...')
% 
% output = sim('simulated_system');
% 
% disp('plotting ...')
% 
% visualize_trajectory(output.tout, output.x, param)
% 
% disp('saving data ...')
% 
% timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
% filename = sprintf('sim_%s.mat', timestamp);
% save(fullfile(pwd, 'simulated_system', 'sim_results', filename), ...
%     'output', 'u', 'param')
% 
% disp('DONE')