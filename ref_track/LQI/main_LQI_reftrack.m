%% LQI reference tracking for the inverted pendulum
% Same controller design as LQI/main_LQI.m.
% The Simulink model (inverted_pendulum_LQI_RT.slx) has been modified so that:
%   - the integrator accumulates (x_ref - x_hat)  instead of (-x_hat)
%   - the gain sees [x_hat-x_ref, v, theta, omega, xi]  instead of [x_hat, ...]
% x_ref_ts (timeseries) is read from the workspace by a From Workspace block.

clear; clc;
h  = 0.01;   % sample time [s]
T  = 15;     % experiment duration [s]
t  = (0:h:T)';
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
run_setup   = 1;   % 1 = build and deploy to hardware, 0 = design only
plot_figures = 1;

proj_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(proj_root, 'global'));

results_dir = fullfile(proj_root, 'ref_track', 'LQI', 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
run_dir = fullfile(results_dir, timestamp);
if ~exist(run_dir, 'dir'), mkdir(run_dir); end

%% 1) Reference trajectory
% Step sequence: 0 -> 0.25 m -> -0.25 m -> 0 m
x_ref_vec = 0.3 * sin(t*1);
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
C_c = [1 0 0 0; 0 0 1 0];
sys_c = ss(A_c, B_c, C_c, zeros(2,1));
sys_d = c2d(sys_c, h);

%% 4) LQI design (same weights as baseline)
C_x   = C_c(1,:);
A_aug = [sys_c.A, zeros(4,1); -C_x, 0];
B_aug = [sys_c.B; 0];
sys_aug_d = c2d(ss(A_aug, B_aug, eye(5), zeros(5,1)), h);

Q_lqi = diag([30, 0.1, 200, 1, 100]);
R_lqi = 10;
K_lqi = dlqr(sys_aug_d.A, sys_aug_d.B, Q_lqi, R_lqi);

fprintf('LQI gain K_lqi:\n'); disp(K_lqi);
fprintf('Closed-loop poles:\n'); disp(eig(sys_aug_d.A - sys_aug_d.B*K_lqi));

%% 5) Kalman filter
R_kf = diag([1e-6, 1e-5]);
Q_kf = diag([1e-8, 1e-3, 1e-12, 1e-2]);
[L, ~] = dlqe(sys_d.A, eye(4), sys_d.C, Q_kf, R_kf);

%% 6) Export to workspace for Simulink
k_pos = p.k_pos; k_neg = p.k_neg;
d_pos = p.d_pos; d_neg = p.d_neg;
Fc = p.Fc;
A  = sys_d.A; B = sys_d.B; C = sys_d.C;

ctrl.h = h; ctrl.T = T;
ctrl.k_pos = k_pos; ctrl.k_neg = k_neg;
ctrl.d_pos = d_pos; ctrl.d_neg = d_neg;
ctrl.Fc = Fc;
ctrl.A = A; ctrl.B = B; ctrl.C = C;
ctrl.K_lqi = K_lqi; ctrl.L = L;
ctrl.Q_lqi = Q_lqi; ctrl.R_lqi = R_lqi;
ctrl.Q_kf  = Q_kf;  ctrl.R_kf  = R_kf;

ts_ctrl = fullfile(run_dir, ['ctrl_' timestamp '.mat']);
save(ts_ctrl, '-struct', 'ctrl');
fprintf('Controller saved to %s\n', ts_ctrl);

%% 7) Run on hardware
if ~run_setup, return; end

mname = 'inverted_pendulum_LQI_RT';
addpath(fullfile(proj_root, 'ref_track', 'LQI', 'inverted-pendulum-inteco'));
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
results_hw.x_hat      = evalin('base', 'x_hat');
results_hw.theta_hat  = evalin('base', 'theta_hat');
results_hw.dx_hat     = evalin('base', 'dx_hat');
results_hw.dtheta_hat = evalin('base', 'dtheta_hat');

res_file = fullfile(run_dir, ['run_hw_' timestamp '.mat']);
save(res_file, '-struct', 'results_hw');
fprintf('Results saved to %s\n', res_file);

%% 9) Plot
if plot_figures
    fig = figure('Name', 'LQI Reference Tracking');

    subplot(3,1,1);
    plot(t, results_hw.x, 'b', t, x_ref_vec, 'r--', 'LineWidth', 1);
    ylabel('x [m]'); title('LQI Reference Tracking'); grid on;
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
    plot_augmented_states(results_hw.x_hat, results_hw.dx_hat, ...
        results_hw.theta_hat, results_hw.dtheta_hat, ...
        h, T, run_dir, timestamp);
end
