function [x_dot, y] = fixed_pendulum_EOM(param, x, u)
% x = [x, v]
% y = sensor data (position)

M   = param.M;
m   = param.m;
c   = param.c;
k_m = param.k_m;

x_dot    = zeros(2,1);
x_dot(1) = x(2);
x_dot(2) = - c/(M+m) * x(2) + k_m/(M+m) * u;

y = x(1);
end