%% LQG + steady-state Kalman design for the inverted pendulum
clear; clc;
h = 0.01;
T = 60;
run_setup = 1;    % 1 = build and run on hardware, 0 = design only
plot_figures = 1; % 1 = save Kalman comparison and augmented state plots after hardware run

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

% add path to the EOM and other functions, assuming this script is in LQG/
proj_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(proj_root, 'global'));

results_dir = fullfile(proj_root, 'LQG', 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
run_dir = fullfile(results_dir, timestamp);
if ~exist(run_dir, 'dir'), mkdir(run_dir); end

%% 1) Load identified parameters
p = load(fullfile(proj_root, 'global', 'param_64_175.mat')).param;
p.g = 9.81;

%% 2) Linearize about the upright equilibrium
% State: z = [x; v; theta; omega], with theta = theta - pi (so theta = 0 is upright)
% Input: F (force on cart, in Newtons). Outputs: y = [x; theta].
% Numerical Jacobian of the continuous EOM with F as input.

f = @(z, F) eom_force(z, F, p.M, p.m, p.b, p.c, p.l);

z_eq = [0; 0; 0; 0];
F_eq = 0;

eps_j = 1e-6;
A = zeros(4,4);
for i = 1:4
    dz = zeros(4,1); dz(i) = eps_j;
    A(:,i) = (f(z_eq+dz, F_eq) - f(z_eq-dz, F_eq)) / (2*eps_j);
end
B = (f(z_eq, F_eq+eps_j) - f(z_eq, F_eq-eps_j)) / (2*eps_j);
C = [1 0 0 0;
     0 0 1 0];
D = zeros(2,1);

sys_c = ss(A, B, C, D);
sys_d = c2d(sys_c, h);
A = sys_d.A; B = sys_d.B; C = sys_d.C; D = sys_d.D;
fprintf('Open-loop poles (continuous):\n'); disp(eig(sys_c.A));
fprintf('Open-loop poles (discrete):\n'); disp(eig(sys_d.A));

%% 3) LQR design (discrete, on augmented system)

Q_lqg = diag([30, 0.1, 200, 1]);   % [x  v  theta  omega]
R_lqg = 3;                         % force penalty


K_lqg = dlqr(sys_d.A, sys_d.B, Q_lqg, R_lqg);   % 1×4

fprintf('LQG gain K_lqg (F = -K_lqg * [x; v; theta; omega]):\n'); disp(K_lqg);
fprintf('Closed-loop poles (discrete):\n');
disp(eig(sys_d.A - sys_d.B * K_lqg));

%% 4) Kalman filter (steady-state, continuous LQE)
% Quantization noise: sigma^2 = delta^2/12
% delta_x = 5.9e-5 m/count, delta_theta = 1.534e-3 rad/count (inferred from data)
R_kf = diag([1e-6, 1e-05]);

% Process noise: treat unmodeled accelerations as white noise entering v
% and omega. Tune q_v, q_w by comparing the estimated velocities against
% a low-pass-filtered finite difference of the measurements.
q_v = 1e-3;   % cart acceleration noise PSD
q_w = 1e-2;   % pendulum angular-acceleration noise PSD
Q_kf = diag([1e-8, q_v, 1e-12, q_w]);

% dlqe: process noise enters all 4 states directly (G = eye(4)).
% Q_kf treated as discrete covariance; scale by h if Q was a continuous PSD.
[L, P_kf] = dlqe(sys_d.A, eye(4), sys_d.C, Q_kf, R_kf);

fprintf('Kalman gain L:\n'); disp(L);
fprintf('Observer poles (A_d - L*C_d):\n'); disp(eig(sys_d.A - L*sys_d.C));

%% 5) Save controller + observer for Simulink 
ctrl.h = h;
ctrl.T = T;
% Dead-zone compensation constants (used in Simulink)
% LQR commands a force F_cmd. Invert the motor map:
%   F_cmd >= 0  ->  u = F_cmd/k_pos + d_pos
%   F_cmd <  0  ->  u = F_cmd/k_neg - d_neg
% A small bias near F_cmd = 0 keeps the cart from chattering; pick a
% threshold like F_dead = 0.02 N if needed.
ctrl.k_pos = p.k_pos; k_pos = p.k_pos;
ctrl.k_neg = p.k_neg; k_neg = p.k_neg;
ctrl.d_pos = p.d_pos; d_pos = p.d_pos;
ctrl.d_neg = p.d_neg; d_neg = p.d_neg;
ctrl.Fc = p.Fc; Fc = p.Fc;

% other stuff to save for Simulink
ctrl.A = sys_d.A; ctrl.B = sys_d.B; ctrl.C = sys_d.C;
ctrl.K_lqg = K_lqg; ctrl.L = L;
ctrl.Q_lqg = Q_lqg; ctrl.R_lqg = R_lqg;
ctrl.Q_kf  = Q_kf;  ctrl.R_kf  = R_kf;
ctrl.sys_c = sys_c; ctrl.sys_d = sys_d; 

out_file = fullfile(results_dir, 'lqr_kalman.mat');
save(out_file, '-struct', 'ctrl');
ts_ctrl_file = fullfile(run_dir, ['ctrl_' timestamp '.mat']);
save(ts_ctrl_file, '-struct', 'ctrl');
fprintf('\nSaved controller + observer to %s\n', out_file);
fprintf('Timestamped copy saved to %s\n', ts_ctrl_file);

%% 6) Run setup or simulate to check the design
if run_setup
    mname = 'inverted_pendulum_LQG';
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
    results_hw.x_hat      = evalin('base', 'x_hat');
    results_hw.theta_hat  = evalin('base', 'theta_hat');
    results_hw.dx_hat     = evalin('base', 'dx_hat');
    results_hw.dtheta_hat = evalin('base', 'dtheta_hat');
    results_hw.F_total    = evalin('base', 'F_total');
    results_hw.delta_F    = evalin('base', 'delta_F');

    res_file = fullfile(run_dir, ['run_hw_' timestamp '.mat']);
    save(res_file, '-struct', 'results_hw');
    fprintf('Hardware results saved to %s\n', res_file);

    if plot_figures
        plot_kalman_comparison(results_hw.x, results_hw.theta, ...
            results_hw.x_hat, results_hw.theta_hat, ...
            results_hw.dx_hat, results_hw.dtheta_hat, ...
            h, T, run_dir, timestamp);
        plot_augmented_states(results_hw.x_hat, results_hw.dx_hat, ...
            results_hw.theta_hat, results_hw.dtheta_hat, ...
            h, T, run_dir, timestamp);
    end
end