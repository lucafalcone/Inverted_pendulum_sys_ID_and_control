function [dx, y] = sysEOM_Coulomb(~, x, u, M, m, b, c, l, k, Fc, varargin)
% x = [x, v, \theta, \omega]
% y = [x, \theta]

g   = 9.81;
epsilon = 30; % Coulomb tanh modifier (to make it continuous)

Delta = M + m*sin(x(3))^2;

num1 = u*k - c*x(2) - Fc*tanh(x(2)*epsilon) - m*l*x(4)^2*sin(x(3)) ...
       - m*g*cos(x(3))*sin(x(3)) - (b/l)*cos(x(3))*x(4);

num2 = m*l*cos(x(3))*num1 - (M+m)*(m*g*l*sin(x(3)) + b*x(4));

dx(1) = x(2);
dx(2) = num1 / Delta;
dx(3) = x(4);
dx(4) = num2 / (m*l^2 * Delta);

y = [x(1), x(3)]';
end