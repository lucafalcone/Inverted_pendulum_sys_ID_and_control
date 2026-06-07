function LTI = build_cw_model(param)

n = param.n;  % Mean motion [rad/s]

%% Continuous-time state-space matrices
% x_dot = Ac * x + Bc * u
LTI.Ac = [0,     0,     1,  0;
          0,     0,     0,  1;
          3*n^2, 0,     0,  2*n;
          0,     0,    -2*n, 0];

LTI.Bc = [0, 0;
          0, 0;
          1, 0;
          0, 1];

LTI.C = eye(param.nx);  % Full state measurement

%% Discretisation via zero-order hold
% x(k+1) = A * x(k) + B * u(k)
Ts = param.Ts;
sys_c = ss(LTI.Ac, LTI.Bc, LTI.C, zeros(param.nx, param.nu));
sys_d = c2d(sys_c, Ts, 'zoh');

LTI.A = sys_d.A;
LTI.B = sys_d.B;

%% Store dimensions
LTI.nx = param.nx;
LTI.nu = param.nu;

end
