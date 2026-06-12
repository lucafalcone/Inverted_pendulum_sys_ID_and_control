%% Swing-up + stabilisation for the inverted pendulum
% Swing-up: Åström-Furuta energy control (Automatica 2000)
%   F = k_e * (E - E0) * sign(ω · cos θ),   E0 = 0 (upright at rest)
% Stabilisation: 4-state discrete LQR once |θ| < 30°
% Kalman filter provides state estimates for both modes.
% θ = 0 is upright throughout; hardware angle origin matches this convention.

clear; clc;
h = 0.01;   % sampling period [s]
T = 30;     % experiment duration [s]  (needs time to swing up)

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

proj_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(proj_root, 'global'));

results_dir = fullfile(proj_root, 'SwingUp', 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
run_dir = fullfile(results_dir, timestamp);
if ~exist(run_dir, 'dir'), mkdir(run_dir); end

%% 1) Physical parameters
p   = load(fullfile(proj_root, 'global', 'param_64_175.mat')).param;
p.g = 9.81;

% Motor / dead-zone constants (also needed by invert_motor block in Simulink)
k_pos = p.k_pos;  k_neg = p.k_neg;
d_pos = p.d_pos;  d_neg = p.d_neg;
Fc    = p.Fc;

% Pendulum physical constants (needed by control_law block in Simulink)
m = p.m;   % mass  [kg]
l = p.l;   % half-length [m]

%% 2) Numerical linearisation about the upright equilibrium (θ = 0)
f   = @(z, F) eom_force(z, F, p.M, p.m, p.b, p.c, p.l);
z0  = [0; 0; 0; 0];
eps = 1e-6;

Ac = zeros(4, 4);
for i = 1:4
    dz = zeros(4,1); dz(i) = eps;
    Ac(:,i) = (f(z0+dz, 0) - f(z0-dz, 0)) / (2*eps);
end
Bc = (f(z0, eps) - f(z0, -eps)) / (2*eps);
Cc = [1 0 0 0;
      0 0 1 0];

sys_c = ss(Ac, Bc, Cc, zeros(2,1));
sys_d = c2d(sys_c, h);
A = sys_d.A;  B = sys_d.B;  C = sys_d.C;

fprintf('Open-loop continuous poles:\n'); disp(eig(Ac).');

%% 3) Discrete LQR  (4-state, no integrator — swing-up does not need cart centering)
% High weight on θ for a fast catch; moderate weight on x to avoid rail hits.
Q_lqr = diag([5, 0.1, 300, 2]);   % [x  v  θ  ω]
R_lqr = 5;
K = dlqr(A, B, Q_lqr, R_lqr);

fprintf('LQR gain K:\n'); disp(K);
fprintf('Closed-loop poles (|z|):\n'); disp(abs(eig(A - B*K)).');

%% 4) Steady-state Kalman filter
% Measurement noise from encoder quantisation
R_kf = diag([1e-6, 1e-5]);            % [σ²_x  σ²_θ]
% Process noise — treat unmodelled acceleration as noise on v and ω
Q_kf = diag([1e-8, 1e-3, 1e-12, 1e-2]);
[L, ~] = dlqe(A, eye(4), C, Q_kf, R_kf);

fprintf('Kalman gain L:\n'); disp(L);
fprintf('Observer poles (|z|):\n'); disp(abs(eig(A - L*C)).');

% Initial Kalman state: pendulum hanging at rest
% (x0_kalman is read by the Kalman DiscreteStateSpace block X0 parameter)
x0_kalman = [0; 0; pi; 0];

%% 5) Energy control parameters
% Åström-Furuta law: F = sat_{Fmax}( k_e · E_err · sign(ω · cos θ) )
% where E_err = ½Jω² + mgl(cos θ − 1),  E_err = 0 at upright at rest.
% Choose k_e so F saturates at Fmax when starting from rest at the bottom
% (E_err = −2mgl at that point).
u_lim = 2.4;                              % leave a 0.1 V margin below safety limit
Fmax  = k_pos * (u_lim - d_pos);         % ≈ 7.3 N
k_e   = Fmax / (2 * m * p.g * l);        % saturates immediately at rest

fprintf('\nEnergy control:  k_e = %.2f N/J,  Fmax = %.2f N\n', k_e, Fmax);
fprintf('(switch to LQR at |θ| < 30°, hardcoded in Simulink control_law block)\n\n');

%% 6) Build and deploy to hardware
% The Simulink model "inverted_pendulum_SU" reads K, L, A, B, C,
% m, l, k_e, Fmax, k_pos, k_neg, d_pos, d_neg, Fc, h, x0_kalman
% directly from the base workspace (this script).
mname = 'inverted_pendulum_SU';
play_run(mname);

fprintf('Waiting for experiment to finish...\n');
while ~strcmp(get_param(mname, 'SimulationStatus'), 'stopped')
    pause(1);
end
fprintf('Experiment complete.\n');

%% 7) Collect hardware results
results_hw.t         = (0:h:T)';
results_hw.x         = evalin('base', 'x');
results_hw.theta     = evalin('base', 'theta');
results_hw.x_hat     = evalin('base', 'x_hat');
results_hw.theta_hat = evalin('base', 'theta_hat');
results_hw.dx_hat    = evalin('base', 'dx_hat');
results_hw.dtheta_hat= evalin('base', 'dtheta_hat');

res_file = fullfile(run_dir, ['run_hw_' timestamp '.mat']);
save(res_file, '-struct', 'results_hw');
fprintf('Hardware results saved to %s\n', res_file);

%% 8) Plots
plot_kalman_comparison(results_hw.x, results_hw.theta, ...
    results_hw.x_hat, results_hw.theta_hat, ...
    results_hw.dx_hat, results_hw.dtheta_hat, ...
    h, T, run_dir, timestamp);
plot_augmented_states(results_hw.x_hat, results_hw.dx_hat, ...
    results_hw.theta_hat, results_hw.dtheta_hat, ...
    h, T, run_dir, timestamp);
