% Parameters
param.M   = 1.0;   % cart mass [kg]
param.m   = 0.2;   % pendulum mass [kg]
param.g   = 9.81;  % gravity [m/s^2]
param.b   = 0.1;   % pendulum damping
param.c   = 0.5;   % cart damping
param.l   = 0.3;   % pendulum length [m]
param.k_m = 1.0;   % motor gain

% Impulse: large force for one short timestep, zero afterwards
dt_imp  = 0.01;          % impulse duration [s]
F_imp   = 10;            % impulse magnitude [N]

u = @(t) F_imp * (t <= dt_imp);

% Initial condition: pendulum hanging DOWN (stable eq.)
% x = [cart_pos, cart_vel, pend_angle, pend_angvel]
x0 = [0; 0; pi; 0];     % pi = hanging down

% Simulate
tspan = [0 5];
[t, x] = ode45(@(t,x) sysEOM(param, x, u(t)), tspan, x0);

% Plot
figure;
subplot(2,1,1)
plot(t, x(:,1)); ylabel('Cart position [m]'); grid on;

subplot(2,1,2)
plot(t, rad2deg(x(:,3))); ylabel('Pendulum angle [deg]'); grid on;
xlabel('Time [s]');
sgtitle('Impulse response');