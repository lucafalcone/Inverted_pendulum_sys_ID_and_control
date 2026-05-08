function [x_dot, y] = fixed_cart_EOM(t, x, u, b_over_m, l, varargin)
% x = [\theta, \omega]
% y = sensor data (position and angle)



g   = 9.81;

x_dot    = zeros(2,1);
x_dot(1) = x(2);
x_dot(2) = -g/l*sin(x(1)) - b_over_m/(l^2) * x(2); % IF theta = 0 around 
% stable equilibrium
% x_dot(2) = g/l*sin(x(1)) - b_over_m/(l^2) * x(2); % IF theta=0 around
% % unstable equilibrium

y = x(1);
end