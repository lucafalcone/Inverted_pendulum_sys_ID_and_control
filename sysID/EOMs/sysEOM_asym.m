function [dx, y] = sysEOM_asym(~, x, u, M, m, b, c, l, k_pos, k_neg, d_pos, d_neg, Fc, varargin)
% x = [x, v, \theta, \omega]
% y = [x, \theta]
% Input map: smooth asymmetric gain + asymmetric dead zone
%   F(u) ~ k_pos*(u - d_pos)  for u >  d_pos
%   F(u) ~ k_neg*(u + d_neg)  for u < -d_neg
%   F(u) ~ 0                  otherwise
% Smoothed via sp(z) = 0.5*(z + sqrt(z^2 + eps_sm^2)), a C^inf max(z,0).
% Cart friction: viscous c*v + Coulomb Fc*tanh(eps_c*v) (smooth sign).

g       = 9.81;
eps_sm  = 1e-3; % smoothing width on the input map  (fixed, not estimated)
eps_c   = 30;   % Coulomb tanh sharpness on cart vel (fixed, not estimated)

sp = @(z) 0.5*(z + sqrt(z.^2 + eps_sm^2));
F  = k_pos * sp(u - d_pos) - k_neg * sp(-(u + d_neg));

Delta = M + m*sin(x(3))^2;

num1 = F - c*x(2) - Fc*tanh(eps_c*x(2)) - m*l*x(4)^2*sin(x(3)) ...
       - m*g*cos(x(3))*sin(x(3)) - (b/l)*cos(x(3))*x(4);

num2 = m*l*cos(x(3))*num1 - (M+m)*(m*g*l*sin(x(3)) + b*x(4));

dx(1) = x(2);
dx(2) = num1 / Delta;
dx(3) = x(4);
dx(4) = num2 / (m*l^2 * Delta);

y = [x(1), x(3)]';
end
