function constraints = constraintgen(LTI, limits, dim )
nx = dim.nx;
nu = dim.nu;
N  = dim.N;

%% 1. Input constraints: -u_max <= u(k) <= u_max for k = 0,...,N-1
% Written as:  [ I; -I ] * U <= [ u_max * ones; u_max * ones ]
A_u = [eye(N*nu); -eye(N*nu)];
b_u = limits.u_max * ones(2*N*nu, 1);
bx0_u = zeros(2*N*nu, nx);

%% 2. State constraints: velocity bounds over prediction horizon
% Bound the approach velocity for safe docking:
%   |dx_dot(k)| <= v_max   and   |dy_dot(k)| <= v_max
% where v_max is derived from the soft-docking requirement.
%
% This constrains predicted states X = T*x0 + S*U.

% Velocity selection matrix: picks velocity states from x
% x = [dx, dy, dx_dot, dy_dot]
F_v = [ 1, 0,  0, 0;     %  dx_dot <= v_max
        -1, 0, 0, 0;     % -dx_dot <= v_max  (i.e. dx_dot >= -v_max)
        0, 0,  0, 0;     %  dy_dot <= v_max
        0, 0,  0, 0];   % -dy_dot <= v_max
e_v = limits.x_max * ones(4, 1);

% Build prediction model
predmod = predmodgen(LTI, dim);

% Stack velocity constraints over the prediction horizon
Fv_bar = kron(eye(N), F_v);
ev_bar = repmat(e_v, N, 1);

A_v   = Fv_bar * predmod.S;
b_v   = ev_bar;
bx0_v = -Fv_bar * predmod.T;

%% Stack all constraints
constraints.A_ineq = [A_u; A_v];
constraints.b_ineq = [b_u; b_v];
constraints.bx0    = [bx0_u; bx0_v];

end
