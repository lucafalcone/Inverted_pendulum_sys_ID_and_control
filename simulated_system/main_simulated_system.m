clear
clc

h = 0.01;
Tsim = 10;
x0 = [0 0 0 0]';

t_u = (0:h:Tsim)';
A = 1;
T_doublet = 1;
u = @(t) A .* (t >= 0 & t < T_doublet) + ...
    -A .* (t >= T_doublet & t < 2*T_doublet);
u = timeseries(u(t_u),t_u);

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