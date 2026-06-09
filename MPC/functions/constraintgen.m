function constraints = constraintgen(LTI, limits, dim)
nx = dim.nx;
nu = dim.nu;
N  = dim.N;

%% 1. Input constraints: -u_max <= u(k) <= u_max for k = 0,...,N-1
A_u   = [eye(N*nu); -eye(N*nu)];
b_u   = limits.u_max * ones(2*N*nu, 1);
bx0_u = zeros(2*N*nu, nx);

%% 2. State constraints: |x(k)| <= x_max + s (soft via single slack s >= 0)
% State vector: [x; v; theta; omega]. First component is cart position.
F_x = [ 1, 0, 0, 0;    %  x <= x_max + s
       -1, 0, 0, 0];   % -x <= x_max + s
e_x = limits.x_max * ones(2, 1);

predmod = predmodgen(LTI, dim);

Fx_bar  = kron(eye(N), F_x);
ex_bar  = repmat(e_x, N, 1);
n_x_rows = 2 * N;

A_x_U  = Fx_bar * predmod.S;    % (2N x N*nu) input part
b_x    = ex_bar;
bx0_x  = -Fx_bar * predmod.T;

%% Stack: augmented variable z = [U (N*nu x 1); s (scalar)]
% Input rows: [A_u,          0] * z <= b_u
% State rows: [A_x_U, -ones  ] * z <= b_x + bx0_x * x0   (slack relaxes bound)
% Slack row:  [0,          -1] * z <= 0                   (enforces s >= 0)
constraints.A_ineq = [A_u,    zeros(2*N*nu, 1);
                      A_x_U, -ones(n_x_rows, 1);
                      zeros(1, N*nu), -1];
constraints.b_ineq = [b_u; b_x; 0];
constraints.bx0    = [bx0_u; bx0_x; zeros(1, nx)];

end
