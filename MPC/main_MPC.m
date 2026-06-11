%% MPC design for the inverted pendulum
% Linearizes about the upright equilibrium, designs the MPC controller with
% terminal LQR cost and discrete Kalman filter, then saves everything for Simulink.
% State: z = [x; v; theta; omega], theta=0 upright. Input: F (force on cart).

clear; clc;
h = 0.05;   % sample time [s]
T = 10;     % experiment duration [s]
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

% MPC tuning
params.Q = diag([50, 0.1, 100, 1]);  % [x  v  theta  omega]
params.R = 5;
limits.u_max = 4;    % max force [N] (~1/3 of the input to the plant)
limits.x_max = 0.2;  % cart rail limit [m]
limits.rho   = 1000; % soft-constraint penalty (increase to enforce harder)
N = 30;              % prediction horizon (steps)

%% Paths
proj_root = fileparts(fileparts(mfilename('fullpath')));
mname = 'inverted_pendulum_MPC';
addpath(fullfile(proj_root, 'MPC', 'functions'));
addpath(fullfile(proj_root, 'global'));
addpath(fullfile(proj_root, 'MPC', 'inverted-pendulum-inteco'));

results_dir = fullfile(proj_root, 'MPC', 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
run_dir = fullfile(results_dir, timestamp);
if ~exist(run_dir, 'dir'), mkdir(run_dir); end

%% 1) Load identified parameters
p = load(fullfile(proj_root, 'global', 'param_64_175.mat')).param;
p.g = 9.81;

%% 2) Linearize about upright equilibrium
f = @(z, F) eom_force(z, F, p.M, p.m, p.b, p.c, p.l);

z_eq = [0; 0; 0; 0];
F_eq = 0;

eps_j = 1e-6;
A_c = zeros(4,4);
for i = 1:4
    dz = zeros(4,1); dz(i) = eps_j;
    A_c(:,i) = (f(z_eq+dz, F_eq) - f(z_eq-dz, F_eq)) / (2*eps_j);
end
B_c = (f(z_eq, F_eq+eps_j) - f(z_eq, F_eq-eps_j)) / (2*eps_j);
C_c = [1 0 0 0;
       0 0 1 0];
D_c = zeros(2,1);

sys_c = ss(A_c, B_c, C_c, D_c);
sys_d = c2d(sys_c, h);
A = sys_d.A; B = sys_d.B; C = sys_d.C;

fprintf('Open-loop poles (continuous):\n'); disp(eig(sys_c.A));
fprintf('Open-loop poles (discrete):\n');   disp(eig(sys_d.A));

%% 3) Terminal cost (LQR solution to DARE)
fprintf('=== Computing terminal ingredients ===\n');
terminal = compute_terminal(sys_d, params);

%% 4) Prediction model and QP cost matrices
weight.Q = params.Q;
weight.R = params.R;
weight.P = terminal.P;

dim.nx = 4;
dim.nu = 1;
dim.N  = N;

predmod    = predmodgen(sys_d, dim);
[H, h_v, ~] = costgen(predmod, weight, dim);
H = blkdiag(H, limits.rho);  % augment for soft-constraint slack variable

%% 5) Constraint matrices
constraints = constraintgen(sys_d, limits, dim);
fprintf('Constraint matrix: A_ineq is %dx%d\n', size(constraints.A_ineq));

%% 6) Discrete steady-state Kalman filter
% Use dlqe (discrete) — same tuning as the working LQI design.
R_kf = diag([1e-6, 1e-05]);
Q_kf = diag([1e-8, 1e-3, 1e-12, 1e-2]);
[L, ~] = dlqe(sys_d.A, eye(4), sys_d.C, Q_kf, R_kf);

fprintf('Kalman gain L:\n'); disp(L);
fprintf('Observer poles (A - L*C):\n'); disp(eig(A - L*C));

%% 7) Save controller + observer for Simulink
% Simulink model reads: A, B, C, L (Kalman block) and H, h_v, constraints,
% dim (MPC MATLAB Function block) plus motor dead-zone constants.
ctrl.h    = h;
ctrl.T    = T;
ctrl.N    = N;
ctrl.k_pos = p.k_pos; k_pos = p.k_pos;
ctrl.k_neg = p.k_neg; k_neg = p.k_neg;
ctrl.d_pos = p.d_pos; d_pos = p.d_pos;
ctrl.d_neg = p.d_neg; d_neg = p.d_neg;
ctrl.Fc = p.Fc; Fc = p.Fc;
ctrl.A    = A; 
ctrl.B    = B;
ctrl.C    = C;
ctrl.L    = L;
ctrl.H    = H;
ctrl.h_v  = h_v;
ctrl.constraints = constraints;
ctrl.dim  = dim;

out_file = fullfile(results_dir, 'mpc_controller.mat');
save(out_file, '-struct', 'ctrl');
ts_ctrl_file = fullfile(run_dir, ['ctrl_' timestamp '.mat']);
save(ts_ctrl_file, '-struct', 'ctrl');
fprintf('\nSaved MPC controller to %s\n', out_file);
fprintf('Timestamped copy saved to %s\n', ts_ctrl_file);

%% GO
play_run(mname);

fprintf('Waiting for experiment to finish...\n');
while ~strcmp(get_param(mname, 'SimulationStatus'), 'stopped')
    pause(1);
end
fprintf('Experiment complete.\n');

% Collect workspace variables written by the To Workspace blocks
results_hw.t          = 0:h:T;
results_hw.x          = evalin('base', 'x');
results_hw.theta      = evalin('base', 'theta');
results_hw.u          = evalin('base', 'outu');
results_hw.x_hat      = evalin('base', 'x_hat');
results_hw.theta_hat  = evalin('base', 'theta_hat');
results_hw.dx_hat     = evalin('base', 'dx_hat');
results_hw.dtheta_hat = evalin('base', 'dtheta_hat');
results_hw.fval       = evalin('base', 'fval');
results_hw.v          = diff(results_hw.x)     ./ h;
results_hw.omega      = diff(results_hw.theta) ./ h;

res_file = fullfile(run_dir, ['run_hw_' timestamp '.mat']);
save(res_file, '-struct', 'results_hw');
fprintf('Hardware results saved to %s\n', res_file);

t_hw   = results_hw.t;
t_diff = t_hw(1:end-1);
fig_hw = figure('Name', 'MPC hardware run');
subplot(5,1,1); plot(t_hw, results_hw.x);
ylabel('x [m]'); title('MPC hardware run'); grid on;
subplot(5,1,2); plot(t_hw, results_hw.theta);
ylabel('\theta [rad]'); grid on;
subplot(5,1,3); plot(t_diff, results_hw.v);
ylabel('v [m/s]'); grid on;
subplot(5,1,4); plot(t_diff, results_hw.omega);
ylabel('\omega [rad/s]'); grid on;
subplot(5,1,5); stairs(t_hw(1:length(results_hw.u)), results_hw.u);
ylabel('F [N]'); xlabel('t [s]'); grid on;
saveas(fig_hw, fullfile(run_dir, ['run_hw_' timestamp '.png']));
fprintf('Hardware plot saved to MPC/results/%s/run_hw_%s.png\n', timestamp, timestamp);

plot_kalman_comparison(results_hw.x, results_hw.theta, ...
    results_hw.x_hat, results_hw.theta_hat, ...
    results_hw.dx_hat, results_hw.dtheta_hat, ...
    h, T, run_dir, timestamp);
