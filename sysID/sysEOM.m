function [x_dot] = sysEOM(param, x, u)

M   = param.M;
m   = param.m;
g   = param.g;
b   = param.b;
c   = param.c;
l   = param.l;
k_m = param.k_m;

f = k_m * u;
D      = (M + m)*m*l^2 - (m*l*cos(x(3)))^2;
phi    = f - c*x(2) + m*l*x(4)^2*sin(x(3));
psi    = m*g*l*sin(x(3)) - b*x(4);

x_dot    = zeros(4,1);
x_dot(1) = x(2);
x_dot(2) = ( m*l^2 * phi  -  m*l*cos(x(3)) * psi ) / D;
x_dot(3) = x(4);
x_dot(4) = ( (M+m) * psi  -  m*l*cos(x(3)) * phi ) / D;

end