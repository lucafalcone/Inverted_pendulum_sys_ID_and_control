clear
clc

h = 0.01;
Tsim = 15;

t_u = (0:h:Tsim)';
u = timeseries(10*sin(2*pi*t_u), t_u);

param.M   = 1.0;   % cart mass [kg]
param.m   = 0.2;   % pendulum mass [kg]
param.g   = 9.81;  % gravity [m/s^2]
param.b   = 0.1;   % pendulum damping
param.c   = 5;   % cart damping
param.l   = 0.3;   % pendulum length [m]
param.k_m = 1.0;   % motor gain

output = sim('simulated_system');

disp('done simulating')

visualize_trajectory(output.tout, output.x, param)

disp('done plotting')