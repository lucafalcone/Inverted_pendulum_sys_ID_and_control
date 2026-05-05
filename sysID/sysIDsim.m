% basically this file is to demonstrate how to achieve parameter
% estimation, in this case it will be error reduction between the sysEOM
% file and the simulink (which uses the same EOM) so ideally the error
% should become 0.

plotting = true;
% define experiment parameters
T = 60;
h = 0.01;
t = (0:h:T)';

% define input signal for the experiment
A = 1;
T_doublet = 1;
u = @(t) A .* (t >= 0 & t < T_doublet) + ...
    -A .* (t >= T_doublet & t < 2*T_doublet);
u = timeseries(u(t),t);

% get real data from real experiment
%TODO
% for now we can get data from a simulation
x0 = [0 0 pi 0]';

param.M   = 0.2;     % cart mass [kg]
param.m   = 0.05;    % pendulum mass [kg]
param.g   = 9.81;    % gravity [m/s^2]
param.b   = 0.001;   % pendulum damping
param.c   = 0.5;   % cart damping
param.l   = 0.2;     % pendulum length [m]
param.k_m = 1.0;     % motor gain

out = simulate_system(h, T, x0, u, param);
y_sim = [out.x(:,1), out.x(:,3)];

% visualize the result
if plotting
    visualize_trajectory(t,out.x,param)
end


