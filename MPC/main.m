clear; clc;
h = 0.01;
T = 10;
run_setup = 1; % set to 0 to run a simulation instead
plot_figures = 1;

% to specify:
params.Q = diag([ 30      0.1     200    1 ]);
params.R = 3;
limits.u_max = 2;
limits.x_max = 1923047;
N = 5; % horizon length

% add path to the EOM and other functions, assuming this script is in LQI/
proj_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(proj_root, 'MPC', 'functions'));

results_dir = fullfile(proj_root, 'MPC', 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

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

%% 3) Compute terminal ingredients
fprintf('=== Computing terminal ingredients ===\n');
terminal = compute_terminal(sys_d, params);

%% 4) Build prediction model and cost function
weight.Q = params.Q;
weight.R = params.R;
weight.P = terminal.P;

dim.nx = 4; % number of states
dim.nu = 1; % number of inputs
dim.N  = N;

predmod = predmodgen(sys_d, dim);
[H, h_v, ~] = costgen(predmod, weight, dim);

%% 5) Build constraint matrices
constraints = constraintgen(sys_d, limits, dim);

%% 4) Kalman filter (steady-state, continuous LQE)
% Measure encoder noise on the stationary rig and put the variances here.
% Until then, start from a sensible guess (encoder resolution ~ 1e-4 m,
% ~ 1e-4 rad standard deviation).
sigma_x     = 1e-4;     % [m]   - REPLACE with measured std of x
sigma_theta = 1e-4;     % [rad] - REPLACE with measured std of theta
R_kf = diag([sigma_x^2, sigma_theta^2]);

% Process noise: treat unmodeled accelerations as white noise entering v
% and omega. Tune q_v, q_w by comparing the estimated velocities against
% a low-pass-filtered finite difference of the measurements.
q_v = 1e-2;   % cart acceleration noise PSD
q_w = 1e-1;   % pendulum angular-acceleration noise PSD
Q_kf = diag([1e-8, q_v, 1e-8, q_w]);

% [kest, L, P] = kalman(...) wants the plant augmented with the noise
% input matrix G. Here process noise drives all states directly: G = I.
G_kf = eye(4);
H_kf = zeros(2,4);
sys_kf = ss(sys_c.A, [sys_c.B G_kf], sys_c.C, [sys_c.D H_kf]);
[~, L, P_kf] = kalman(sys_kf, Q_kf, R_kf);

fprintf('Kalman gain L:\n'); disp(L);
fprintf('Observer poles (A - L*C):\n'); disp(eig(A - L*C));



%% GO
mname = 'inverted_pendulum_MPC';
play_run(mname)