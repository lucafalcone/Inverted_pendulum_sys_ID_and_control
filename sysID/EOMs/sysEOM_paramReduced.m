function [x_dot, y] = sysEOM_paramReduced(~, x, u, M, m, c, k_m, varargin)
% x = [x, v, \theta, \omega]
% y = [x, \theta]

l = 0.379670764674439;
b_over_m = 0.005966407584023;

g   = 9.81;
f   = -k_m * u;

D   = (M + m)*m*l^2 - (m*l*cos(x(3)))^2;
phi = f - c*x(2) - m*l*x(4)^2*sin(x(3));       
psi = m*g*l*sin(x(3)) + b_over_m * m *x(4);               

x_dot    = zeros(4,1);
x_dot(1) = x(2);
x_dot(2) = (  m*l^2  * phi  -  m*l*cos(x(3)) * psi ) / D;
x_dot(3) = x(4);
x_dot(4) = -( (M+m)  * psi  -  m*l*cos(x(3)) * phi ) / D;  

y = [x(1), x(3)]';
end