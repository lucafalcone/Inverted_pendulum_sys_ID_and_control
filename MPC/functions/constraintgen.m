function constraints = constraintgen(LTI, limits, dim)
nx = dim.nx;
nu = dim.nu;
N  = dim.N;

%% 1. Input constraints: -u_max <= u(k) <= u_max for k = 0,...,N-1
A_u   = [eye(N*nu); -eye(N*nu)];
b_u   = limits.u_max * ones(2*N*nu, 1);
bx0_u = zeros(2*N*nu, nx);

%% 2. State constraints: cart position |x(k)| <= x_max over horizon
% State vector: [x; v; theta; omega]. First component is cart position.
F_x = [ 1, 0, 0, 0;    %  x <= x_max
       -1, 0, 0, 0];   % -x <= x_max
e_x = limits.x_max * ones(2, 1);

predmod = predmodgen(LTI, dim);

Fx_bar  = kron(eye(N), F_x);
ex_bar  = repmat(e_x, N, 1);

A_x   = Fx_bar * predmod.S;
b_x   = ex_bar;
bx0_x = -Fx_bar * predmod.T;

%% Stack all constraints
constraints.A_ineq = [A_u; A_x];
constraints.b_ineq = [b_u; b_x];
constraints.bx0    = [bx0_u; bx0_x];

end
