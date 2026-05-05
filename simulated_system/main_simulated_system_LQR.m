clear
clc

h = 0.01;
Tsim = 20;
x0 = [0 0.1 0.2 -0.1]';

t_u = (0:h:Tsim)';
v = timeseries(2*sin(2*pi*t_u) + 1*sin(10*pi*t_u), t_u);

param.M   = 0.2;     % cart mass [kg]
param.m   = 0.05;    % pendulum mass [kg]
param.g   = 9.81;    % gravity [m/s^2]
param.b   = 0.001;   % pendulum damping
param.c   = 0.005;   % cart damping
param.l   = 0.2;     % pendulum length [m]
param.k_m = 1.0;     % motor gain

Q = diag([10, 1, 100, 1]);
R = 1;
load('linearized_system')
K = -lqr(linsys1, Q, R);

disp('simulating ...')

output = sim('simulated_system_LQR');

disp('plotting ...')

visualize_trajectory(output.tout, output.x, param)

disp('saving data ...')

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
filename = sprintf('simLQR_%s.mat', timestamp);
save(fullfile(pwd, 'simulated_system', 'sim_results', filename), ...
    'output', 'v', 'param', 'Q', 'R', 'K')

disp('DONE')