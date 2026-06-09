%% Black-box closed-loop system identification of the inverted pendulum
% Identifies a 4th-order discrete linear state-space model from closed-loop
% hardware data using OKID+ERA (no equations of motion are used).
%
% Workflow:
%   1) Load bootstrap LQR (from LQI/main.m) so the hardware can run.
%   2) Run Simulink experiment: LQR stabilises the pendulum upright while a
%      chirp reference perturbation excites the closed-loop dynamics.
%      Simulink must log: x, theta, outu, t  to the workspace.
%   3) Load closed-loop I/O data {outu(t), [x(t), theta(t)]}.
%   4) OKID+ERA: ARX-style closed-loop regression → system Markov parameters
%      → block-Hankel SVD → minimal state-space realisation A,B,C.
%   5) LQR + Kalman redesign on the identified model.
%   6) Observer-based prediction validation (K-step-ahead, K=10).
%   7) Bode / pole-zero comparison vs physics-based linearisation.
%   8) Save bb_identified.mat for Simulink.

clear; clc;
h  = 0.01;
T  = 10;
run_setup    = 1;   % 1: run hardware; 0: skip to identification
plot_figures = 1;
nx = 4;             % model order: [x, v, theta, omega]
p  = 100;           % Markov parameters for OKID (covers 1 s of impulse response)
K_hor = 10;         % prediction horizon for validation [samples]

t = 0:h:T;

proj_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(proj_root, 'black_box', 'functions'));
addpath(fullfile(proj_root, 'global'));

results_dir = fullfile(proj_root, 'black_box', 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% 1) Bootstrap: load the stabilising LQR so the hardware experiment is possible
ctrl_file = fullfile(results_dir, 'lqr_kalman.mat');
if ~exist(ctrl_file, 'file')
    error('Bootstrap controller not found. Run LQI/main.m first.');
end
ctrl = load(ctrl_file);
A = ctrl.A; B = ctrl.B; C = ctrl.C;
K_lqi = ctrl.K_lqi; L = ctrl.L;
k_pos = ctrl.k_pos; k_neg = ctrl.k_neg;
d_pos = ctrl.d_pos; d_neg = ctrl.d_neg;
fprintf('Bootstrap LQR loaded.\n');
fprintf('  Open-loop poles:   '); fprintf('% .4f  ', eig(ctrl.A)'); fprintf('\n');
fprintf('  Closed-loop poles: '); fprintf('% .4f  ', eig(ctrl.A - ctrl.B*ctrl.K_lqi(1:4))'); fprintf('\n\n');

%% 2) Hardware experiment
% The Simulink model runs the LQR while a chirp (e.g. 0.1–5 Hz, ±0.02 m)
% is added to the cart-position reference to persistently excite the system.
% After the run, Simulink must have logged x, theta, outu, t to the workspace.
if run_setup
    mname = 'inverted_pendulum_LQRBB';
    play_run(mname);
    pause(T + 2);
    save(fullfile(results_dir, 'cl_experiment.mat'), 'x', 'theta', 'outu', 't');
end

%% 3) Load closed-loop data
% Priority: saved hardware data; fallback: sysID chirp files collected
% with the LQR running (these have theta in ±0.5 rad around the upright).
exp_file = fullfile(results_dir, 'cl_experiment.mat');
if exist(exp_file, 'file')
    E = load(exp_file);
    x_cl = E.x(:); th_cl = E.theta(:); u_cl = E.outu(:);
    fprintf('Hardware data: %d samples (%.1f s)\n', numel(x_cl), numel(x_cl)*h);
else
    fallback = {
        fullfile(proj_root, 'sysID', 'real_data', 'chirp2_20260608_113815.mat')
        fullfile(proj_root, 'sysID', 'real_data', 'chirp2_20260608_113905.mat')
    };
    x_cl = []; th_cl = []; u_cl = [];
    for i = 1:numel(fallback)
        S     = load(fallback{i});
        x_cl  = [x_cl;  S.x(:)];    %#ok<AGROW>
        th_cl = [th_cl; S.theta(:)]; %#ok<AGROW>
        u_cl  = [u_cl;  S.outu(:)];  %#ok<AGROW>
    end
    fprintf('Fallback sysID data: %d samples, theta in [%.3f, %.3f] rad\n', ...
        numel(x_cl), min(th_cl), max(th_cl));
end
N_samp = numel(x_cl);

%% 4) OKID+ERA identification
% Remove DC: identification is of the perturbation dynamics around the
% operating point, not absolute position.
y_id = [x_cl - mean(x_cl), th_cl - mean(th_cl)];
u_id = u_cl - mean(u_cl);

n_tr = floor(0.7 * N_samp);
u_tr = u_id(1:n_tr);     y_tr = y_id(1:n_tr, :);
u_vl = u_id(n_tr+1:end); y_vl = y_id(n_tr+1:end, :);

fprintf('OKID+ERA: nx=%d, p=%d, %d training / %d validation samples\n', ...
    nx, p, n_tr, N_samp - n_tr);
[A_bb, B_bb, C_bb, sv] = okid_era(u_tr, y_tr, nx, p);
sys_bb = ss(A_bb, B_bb, C_bb, zeros(2,1), h);

fprintf('\nIdentified open-loop poles:\n');
disp(eig(A_bb));
if any(abs(eig(A_bb)) > 1)
    fprintf('  -> Unstable pole confirmed (expected for upright pendulum)\n\n');
end

if plot_figures
    fig_sv = figure('Name', 'Hankel singular values');
    semilogy(1:min(15, length(sv)), sv(1:min(15, length(sv))), 'o-');
    xlabel('Index'); ylabel('Singular value (log scale)');
    title('Hankel singular values — gap confirms model order');
    xline(nx, 'r--', sprintf('nx = %d', nx), 'LabelHorizontalAlignment', 'left');
    grid on;
    saveas(fig_sv, fullfile(results_dir, 'hankel_sv.png'));
end

%% 5) LQR + Kalman design on identified model
% LQR weights match the bootstrap controller so performance is comparable.
Q_lqr = diag([30, 0.1, 200, 1]);
R_lqr = 1000;
try
    K_bb = dlqr(A_bb, B_bb, Q_lqr, R_lqr);
    fprintf('LQR redesigned on identified model.\n');
    fprintf('  CL poles: '); fprintf('% .4f  ', eig(A_bb - B_bb*K_bb)'); fprintf('\n');
catch e_lqr
    warning('dlqr on identified model failed (%s). Using bootstrap gains.');
    K_bb = ctrl.K_lqi(1:4);
end

R_kf = diag([1e-6, 1e-5]);
Q_kf = diag([1e-8, 1e-3, 1e-12, 1e-2]);
try
    [L_bb, ~] = dlqe(A_bb, eye(4), C_bb, Q_kf, R_kf);
    fprintf('Kalman redesigned on identified model.\n');
    fprintf('  Observer poles: '); fprintf('% .4f  ', eig(A_bb - L_bb*C_bb)'); fprintf('\n\n');
catch e_kf
    warning('dlqe on identified model failed (%s). Using bootstrap observer.');
    L_bb = ctrl.L;
end

%% 6) Validation: K-step-ahead observer prediction
% Run the identified model as a Luenberger observer, then simulate K_hor
% steps open-loop before re-injecting the measurement.  For an unstable
% system this is the only meaningful finite-horizon metric.
N_vl = size(u_vl, 1);
y_pred = zeros(N_vl, 2);
x_obs  = zeros(nx, 1);

for k = 1:N_vl
    if mod(k-1, K_hor) == 0
        % Re-inject measurement: correct state via observer update
        y_k = y_vl(k, :)';
        x_obs = x_obs + L_bb * (y_k - C_bb * x_obs);
    end
    y_pred(k, :) = (C_bb * x_obs)';
    x_obs = A_bb * x_obs + B_bb * u_vl(k);
end

fit_vl = zeros(1, 2);
for i = 1:2
    e = y_vl(:,i) - y_pred(:,i);
    fit_vl(i) = 100 * max(0, 1 - norm(e) / norm(y_vl(:,i) - mean(y_vl(:,i))));
end
fprintf('%d-step-ahead prediction fit (validation):\n', K_hor);
fprintf('  x: %.1f%%   theta: %.1f%%\n\n', fit_vl(1), fit_vl(2));

if plot_figures
    t_vl = (0:N_vl-1)' * h;
    fig_val = figure('Name', 'Black-box validation');
    subplot(2,1,1);
    plot(t_vl, y_vl(:,1), 'b', t_vl, y_pred(:,1), 'r--');
    ylabel('x [m]'); grid on;
    title(sprintf('Identified model: %d-step-ahead prediction vs measured (validation set)', K_hor));
    legend('measured', sprintf('model (fit %.1f%%)', fit_vl(1)));
    subplot(2,1,2);
    plot(t_vl, y_vl(:,2), 'b', t_vl, y_pred(:,2), 'r--');
    ylabel('\theta [rad]'); xlabel('t [s]'); grid on;
    legend('measured', sprintf('model (fit %.1f%%)', fit_vl(2)));
    saveas(fig_val, fullfile(results_dir, 'bb_validation.png'));
end

%% 7) Bode and pole-zero comparison vs physics-based linearisation
% Scale B_phys by k_pos so the physics model input is in motor-command units,
% matching the black-box model (which also has outu as input).
p_phys = load(fullfile(proj_root, 'global', 'param_64_175.mat')).param;
p_phys.g = 9.81;
f_eom = @(z, F) eom_force(z, F, p_phys.M, p_phys.m, p_phys.b, p_phys.c, p_phys.l);
eps_j = 1e-6; z_eq = zeros(4,1);
A_ph = zeros(4,4);
for i = 1:4
    dz = zeros(4,1); dz(i) = eps_j;
    A_ph(:,i) = (f_eom(z_eq+dz,0) - f_eom(z_eq-dz,0)) / (2*eps_j);
end
B_ph = (f_eom(z_eq,eps_j) - f_eom(z_eq,-eps_j)) / (2*eps_j);
sys_phys_d = c2d(ss(A_ph, B_ph * p_phys.k_pos, [1 0 0 0; 0 0 1 0], zeros(2,1)), h);

fprintf('Physics model poles (discrete):\n');
disp(eig(sys_phys_d.A));

if plot_figures
    w = logspace(-1, 2.5, 400);
    fig_bode = figure('Name', 'Bode comparison');
    bodeplot(sys_bb, sys_phys_d, w);
    legend('OKID+ERA (black-box)', 'Physics linearisation', 'Location', 'best');
    title('Bode: black-box identified vs physics-based linearisation');
    saveas(fig_bode, fullfile(results_dir, 'bode_comparison.png'));

    fig_pz = figure('Name', 'Pole-zero comparison');
    pzmap(sys_bb, sys_phys_d);
    legend('OKID+ERA (black-box)', 'Physics linearisation');
    title('Pole-zero map: identified vs physics model');
    grid on;
    saveas(fig_pz, fullfile(results_dir, 'pzmap_comparison.png'));
end

%% 8) Save results
out.h       = h;        out.T       = T;
out.A       = A_bb;     out.B       = B_bb;     out.C       = C_bb;
out.K_lqi   = K_bb;     out.L       = L_bb;
out.Q_lqr   = Q_lqr;    out.R_lqr   = R_lqr;
out.Q_kf    = Q_kf;     out.R_kf    = R_kf;
out.sv      = sv;
out.fit_val = fit_vl;
out.k_pos   = p_phys.k_pos;   out.k_neg = p_phys.k_neg;
out.d_pos   = p_phys.d_pos;   out.d_neg = p_phys.d_neg;
out.sys_bb = ss(A_bb, B_bb, C_bb, zeros(2,1)); 

save(fullfile(results_dir, 'bb_identified.mat'), '-struct', 'out');
fprintf('Saved: black-box model + controller → bb_identified.mat\n');

%% 9) Kalman comparison plot (after hardware run)
if plot_figures && run_setup && exist('x_hat', 'var')
    plot_kalman_comparison(x, theta, x_hat, theta_hat, dx_hat, dtheta_hat, h, T, results_dir);
end
