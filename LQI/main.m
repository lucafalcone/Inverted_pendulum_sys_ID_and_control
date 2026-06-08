%% LQI + steady-state Kalman design for the inverted pendulum
% Linearizes sysEOM_asym about the UPRIGHT equilibrium,
% with the cart FORCE F as the control input (motor dead zone inverted
% downstream in Simulink). Exports K (1×5 LQI gain), L, A, B, C and the
% dead-zone compensation constants for the Simulink model.
% The 5th state is xi = integral(x_ref - x_hat), cart position error.

clear; clc;
h = 0.01;
T = 30;
run_setup = 1; % set to 0 to run a simulation instead

% add path to the EOM and other functions, assuming this script is in LQI/
proj_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(proj_root, 'LQI', 'functions'));

%% 1) Load identified parameters
p = load(fullfile(proj_root, 'sysID', 'results', 'param_60_829.mat')).param;
p.g = 9.81;

%% 2) Linearize about the upright equilibrium
% State: z = [x; v; phi; omega], with phi = theta - pi (so phi = 0 is upright)
% Input: F (force on cart, in Newtons). Outputs: y = [x; phi].
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
fprintf('Open-loop poles (continuous):\n'); disp(eig(sys_c.A));
fprintf('Open-loop poles (discrete):\n'); disp(eig(sys_d.A));

%% 3) LQI design (discrete, on augmented system)
% Augment the 4-state plant with xi = integral(x_ref - x_hat).
% xi_dot = x_ref - x = -x  (x_ref = 0), so the integrator removes
% steady-state cart-position error without changing the pendulum design.
C_x   = C(1,:);                    % [1 0 0 0] — selects cart position
A_aug = [A,    zeros(4,1);
         -C_x, 0          ];       % 5×5 augmented A
B_aug = [B; 0];                    % 5×1 augmented B

sys_aug_c = ss(A_aug, B_aug, eye(5), zeros(5,1));
sys_aug_d = c2d(sys_aug_c, h);     % ZOH discretisation

% Tune weights. The last diagonal entry weights the integral state xi.
Q_lqi = diag([30, 0.1, 200, 1, 5]);   % [x  v  theta  omega  xi]
R_lqi = 3;                             % force penalty (same as before)

K_lqi = dlqr(sys_aug_d.A, sys_aug_d.B, Q_lqi, R_lqi);   % 1×5

fprintf('LQI gain K_lqi (F = -K_lqi * [x; v; theta; omega; xi]):\n'); disp(K_lqi);
fprintf('Closed-loop poles (augmented discrete):\n');
disp(eig(sys_aug_d.A - sys_aug_d.B * K_lqi));

%% 4) Kalman filter (steady-state, continuous LQE)
% Quantization noise: sigma^2 = delta^2/12
% delta_x = 5.9e-5 m/count, delta_theta = 1.534e-3 rad/count (inferred from data)
R_kf = diag([2.95e-10, 1.96e-7]);

% Process noise: treat unmodeled accelerations as white noise entering v
% and omega. Tune q_v, q_w by comparing the estimated velocities against
% a low-pass-filtered finite difference of the measurements.
q_v = 1e-2;   % cart acceleration noise PSD
q_w = 1e-1;   % pendulum angular-acceleration noise PSD
Q_kf = diag([1e-8, q_v, 1e-8, q_w]);

% [kest, L, P] = kalman(...) wants the plant augmented with the noise
% input matrix G. Here process noise drives all states directly: G = I.
G = eye(4);
H = zeros(2,4);
sys_kf = ss(A, [B G], C, [D H]);
[~, L, P_kf] = kalman(sys_kf, Q_kf, R_kf);

fprintf('Kalman gain L:\n'); disp(L);
fprintf('Observer poles (A - L*C):\n'); disp(eig(A - L*C));

%% 5) Save controller + observer for Simulink 
ctrl.h = h;
ctrl.T = T;
% Dead-zone compensation constants (used in Simulink)
% LQR commands a force F_cmd. Invert the motor map:
%   F_cmd >= 0  ->  u = F_cmd/k_pos + d_pos
%   F_cmd <  0  ->  u = F_cmd/k_neg - d_neg
% A small bias near F_cmd = 0 keeps the cart from chattering; pick a
% threshold like F_dead = 0.02 N if needed.
ctrl.k_pos = p.k_pos;
ctrl.k_neg = p.k_neg;
ctrl.d_pos = p.d_pos;
ctrl.d_neg = p.d_neg;

% other stuff to save for Simulink
ctrl.A = A; ctrl.B = B; ctrl.C = C;
ctrl.K = K_lqi; ctrl.L = L;
ctrl.Q_lqi = Q_lqi; ctrl.R_lqi = R_lqi;
ctrl.Q_kf  = Q_kf;  ctrl.R_kf  = R_kf;

out_file = fullfile(proj_root, 'LQI', 'lqr_kalman.mat');
save(out_file, '-struct', 'ctrl');
fprintf('\nSaved controller + observer to %s\n', out_file);

%% 6) Run setup or simulate to check the design
if run_setup
    mname = 'inverted_pendulum_template';
    play_run(mname)
else
    % Quick sanity sim of the linear closed loop (discrete LQI + Kalman)
    t  = 0:h:5;
    N  = length(t);
    z0 = [0; 0; 0.15; 0];          % 0.15 rad ~ 8.6 deg initial tilt

    Ad = sys_d.A; Bd = sys_d.B; Cd = sys_d.C;
    L_d = dlqe(Ad, eye(4), Cd, Q_kf, R_kf);

    K1 = K_lqi(1:4);   % gain on the 4 plant states
    K5 = K_lqi(5);     % integral gain

    q     = zeros(4, N);  q(:,1) = z0;   % true plant state
    q_hat = zeros(4, N);                  % Kalman estimate (starts at zero = full error)
    xi    = zeros(1, N);                  % cart-position integrator state

    for k = 1:N-1
        y_k         = Cd * q(:,k);                            % noiseless measurement
        F_k         = -K1 * q_hat(:,k) - K5 * xi(k);         % LQI control law
        q(:,k+1)    = Ad * q(:,k) + Bd * F_k;                 % true dynamics
        q_hat(:,k+1)= (Ad - L_d*Cd)*q_hat(:,k) + Bd*F_k + L_d*y_k;  % Kalman
        xi(k+1)     = xi(k) + h * (-q_hat(1,k));              % xi_dot = x_ref - x_hat
    end

    figure('Name','Linear LQI+KF closed loop');
    subplot(3,1,1); plot(t, q(1,:)); ylabel('x [m]'); grid on;
    title('Linear closed-loop response from initial tilt (LQI + Kalman)');
    subplot(3,1,2); plot(t, q(3,:)); ylabel('\theta [rad]'); grid on;
    subplot(3,1,3); plot(t, xi);     ylabel('\xi [m·s]'); xlabel('t [s]'); grid on;
    title('Integrator state \xi');
end