% the following does not account for any Coulomb
% friction/dead-zone/asymmetry. It is used just as analytical
% validation for the numerically obtained matrices.

M = param.M;
m = param.m;
l = param.l;
c = param.c;
b = param.b;
k = (param.k_pos + param.k_neg)/2;
g = 9.81;

A= [0,       1,             0,                0;
    0,    -c/M,        -g*m/M,          b/(M*l);
    0,       0,             0,                1;
    0, c/(M*l), g*(M+m)/(M*l), -b*(M+m)/(M*m*l)];
B = [0; k/M; 0; -k/(M*l)];
C = [1 0 0 0;
    0 0 1 0];
D = 0;

% the matrices are physically consistent, they differ mostly in the terms
% a22 and a42, because of the presence of the parameter c.
% most probably the dead-zones and Coulomb appear in the matrices in the
% same place as the parameter c, since they have similar effects on the
% cart (opposing to the input);