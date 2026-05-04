function [x_dot, y] = fixed_cart_EOM(param, x)
% x = [\theta, \omega]
% y = sensor data (position and angle)

m   = param.m;
g   = param.g;
b   = param.b;
l   = param.l;

x_dot    = zeros(2,1);
x_dot(1) = x(2);
x_dot(2) = g/l*sin(x(1)) - b/m/(l^2) * x(2);

y = x(1);
end