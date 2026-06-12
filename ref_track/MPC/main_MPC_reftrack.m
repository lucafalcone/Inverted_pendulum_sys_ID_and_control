%% MPC reference tracking for the inverted pendulum
% Same QP design as MPC/main_MPC.m.
% The Simulink model (inverted_pendulum_MPC_RT.slx) passes the shifted state
%   x0_shifted = x_hat - [x_ref; 0; 0; 0]
% to the MPC MATLAB Function block. Because [x_ref; 0; 0; 0] is a fixed
% point of A_d (A_d * [x_ref;0;0;0] = [x_ref;0;0;0]), this shift is exactly
% equivalent to the full tracking QP formulation with no approximation.
% x_ref_ts (timeseries) is read from the workspace by a From Workspace block.

clear; clc;
h  = 0.03;   % sample time [s]
T  = 15;     % experiment duration [s]
t  = (0:h:T)';
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

proj_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
mname = 'inverted_pendulum_MPC_RT';
addpath(fullfile(proj_root, 'MPC', 'functions'));
addpath(fullfile(proj_root, 'global'));
addpath(fullfile(proj_root, 'ref_track', 'MPC', 'inverted-pendulum-inteco'));

results_dir = fullfile(proj_root, 'ref_track', 'MPC', 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
run_dir = fullfile(results_dir, timestamp);
if ~exist(run_dir, 'dir'), mkdir(run_dir); end

%% 1) Reference trajectory
% Step sequence: 0 -> 0.25 m -> -0.25 m -> 0 m
x_ref_vec = zeros(size(t));
% x_ref_vec(t >= 3  & t < 7 ) =  0.25;
% x_ref_vec(t >= 7  & t < 11) = -0.25;
% x_ref_vec(t >= 11)          =  0;
% x_ref_ts = timeseries(x_ref_vec, t);
x_ref_vec = 0.3 * sin(t*0.6);
x_ref_ts = timeseries(x_ref_vec, t);

%% 2) Load identified parameters
p = load(fullfile(proj_root, 'global', 'param_64_175.mat')).param;
p.g = 9.81;

%% 3) Linearise about upright equilibrium
f = @(z, F) eom_force(z, F, p.M, p.m, p.b, p.c, p.l);
z_eq = [0;0;0;0]; F_eq = 0; eps_j = 1e-6;
A_c = zeros(4,4);
for i = 1:4
    dz = zeros(4,1); dz(i) = eps_j;
    A_c(:,i) = (f(z_eq+dz, F_eq) - f(z_eq-dz, F_eq)) / (2*eps_j);
end
B_c = (f(z_eq, F_eq+eps_j) - f(z_eq, F_eq-eps_j)) / (2*eps_j);
sys_c = ss(A_c, B_c, [1 0 0 0; 0 0 1 0], zeros(2,1));
sys_d = c2d(sys_c, h);
A = sys_d.A; B = sys_d.B; C = sys_d.C;

fprintf('Open-loop poles (continuous):\n'); disp(eig(sys_c.A));

%% 4) MPC design (same as baseline)
params.Q = diag([50, 0.1, 100, 1]);
params.R = 5;
limits.u_max = 4;
limits.x_max = 0.5;
limits.rho   = 1000;
N = 30;

dim.nx = 4; dim.nu = 1; dim.N = N;

terminal = compute_terminal(sys_d, params);
weight.Q = params.Q; weight.R = params.R; weight.P = terminal.P;

predmod    = predmodgen(sys_d, dim);
[H, h_v, ~] = costgen(predmod, weight, dim);
H = blkdiag(H, limits.rho);

constraints = constraintgen(sys_d, limits, dim);
fprintf('Constraint matrix: A_ineq is %dx%d\n', size(constraints.A_ineq));

%% 5) Kalman filter
R_kf = diag([1e-6, 1e-5]);
Q_kf = diag([1e-8, 1e-3, 1e-12, 1e-2]);
[L, ~] = dlqe(sys_d.A, eye(4), sys_d.C, Q_kf, R_kf);
fprintf('Kalman gain L:\n'); disp(L);

%% 6) Export to workspace for Simulink
k_pos = p.k_pos; k_neg = p.k_neg;
d_pos = p.d_pos; d_neg = p.d_neg;
Fc = p.Fc;

ctrl.h = h; ctrl.T = T; ctrl.N = N;
ctrl.k_pos = k_pos; ctrl.k_neg = k_neg;
ctrl.d_pos = d_pos; ctrl.d_neg = d_neg;
ctrl.Fc = Fc;
ctrl.A = A; ctrl.B = B; ctrl.C = C;
ctrl.L = L; ctrl.H = H; ctrl.h_v = h_v;
ctrl.constraints = constraints; ctrl.dim = dim;

ts_ctrl = fullfile(run_dir, ['ctrl_' timestamp '.mat']);
save(ts_ctrl, '-struct', 'ctrl');
fprintf('Controller saved to %s\n', ts_ctrl);

%% 7) Run on hardware
play_run(mname);

fprintf('Waiting for experiment to finish...\n');
while ~strcmp(get_param(mname, 'SimulationStatus'), 'stopped')
    pause(1);
end
fprintf('Experiment complete.\n');

%% 8) Collect and save results
results_hw.t          = t;
results_hw.x_ref      = x_ref_vec;
results_hw.x          = evalin('base', 'x');
results_hw.theta      = evalin('base', 'theta');
results_hw.u          = evalin('base', 'outu');
results_hw.x_hat      = evalin('base', 'x_hat');
results_hw.theta_hat  = evalin('base', 'theta_hat');
results_hw.dx_hat     = evalin('base', 'dx_hat');
results_hw.dtheta_hat = evalin('base', 'dtheta_hat');
results_hw.fval       = evalin('base', 'fval');

res_file = fullfile(run_dir, ['run_hw_' timestamp '.mat']);
save(res_file, '-struct', 'results_hw');
fprintf('Results saved to %s\n', res_file);

%% 9) Plot
fig = figure('Name', 'MPC Reference Tracking');

subplot(3,1,1);
plot(t, results_hw.x, 'b', t, x_ref_vec, 'r--', 'LineWidth', 1);
ylabel('x [m]'); title('MPC Reference Tracking'); grid on;
legend('x_{hat}', 'x_{ref}'); xlim([0 T]);

subplot(3,1,2);
plot(t, results_hw.theta_hat, 'b', 'LineWidth', 1);
ylabel('\theta [rad]'); grid on; xlim([0 T]);

subplot(3,1,3);
plot(t, results_hw.x_hat - x_ref_vec, 'k', 'LineWidth', 1);
ylabel('tracking error [m]'); xlabel('t [s]'); grid on; xlim([0 T]);

saveas(fig, fullfile(run_dir, ['reftrack_' timestamp '.png']));

plot_kalman_comparison(results_hw.x, results_hw.theta, ...
    results_hw.x_hat, results_hw.theta_hat, ...
    results_hw.dx_hat, results_hw.dtheta_hat, ...
    h, T, run_dir, timestamp);
