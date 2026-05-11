%% LQR + steady-state Kalman design for the inverted pendulum
% Linearizes sysEOM_asym about the UPRIGHT equilibrium (theta = pi),
% with the cart FORCE F as the control input (motor dead zone inverted
% downstream in Simulink). Exports K, L, A, B, C and the dead-zone
% compensation constants for the Simulink model.

clear; clc;
h = 0.01;
T = 30;

%% 1) Load identified parameters
S = load(fullfile('Results', 'param_60_829.mat'));
% file stores either a struct 'param' or the raw fields - handle both
if isfield(S, 'param')
    p = S.param;
else
    p = S;
end

M     = p.M;
m     = p.m;
b     = p.b;     % pendulum pivot viscous friction
c     = p.c;     % cart viscous friction
l     = p.l;
k_pos = p.k_pos;
k_neg = p.k_neg;
d_pos = p.d_pos;
d_neg = p.d_neg;
Fc    = p.Fc;    % Coulomb cart friction (ignored in linearization)

g = 9.81;

%% 2) Linearize about the upright equilibrium
% State: z = [x; v; phi; omega], with phi = theta - pi (so phi = 0 is upright)
% Input: F (force on cart, in Newtons). Outputs: y = [x; phi].
% Numerical Jacobian of the continuous EOM with F as input.

f = @(z, F) eom_force(z, F, M, m, b, c, l);

z_eq = [0; 0; 0; 0];   % upright, at rest, at origin (in phi coordinates)
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
fprintf('Open-loop poles (continuous):\n'); disp(eig(A));

%% 3) LQR design (continuous)
% Tune these. Heavy weight on phi keeps the pendulum upright; small weight
% on x prevents drift but slowly. R penalizes force effort.
Q_lqr = diag([ 30      0.1     200     1   ]);   % [x  v  phi  omega]
R_lqr =  3;                                    % force penalty

K = lqr(A, B, Q_lqr, R_lqr);
fprintf('LQR gain K (force = -K * [x;v;phi;omega]):\n'); disp(K);
fprintf('Closed-loop poles (A - B*K):\n'); disp(eig(A - B*K));

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
G = eye(4);
H = zeros(2,4);
sys_kf = ss(A, [B G], C, [D H]);
[~, L, P_kf] = kalman(sys_kf, Q_kf, R_kf);

fprintf('Kalman gain L:\n'); disp(L);
fprintf('Observer poles (A - L*C):\n'); disp(eig(A - L*C));

%% 5) Dead-zone compensation constants (used in Simulink)
% LQR commands a force F_cmd. Invert the motor map:
%   F_cmd >= 0  ->  u = F_cmd/k_pos + d_pos
%   F_cmd <  0  ->  u = F_cmd/k_neg - d_neg
% A small bias near F_cmd = 0 keeps the cart from chattering; pick a
% threshold like F_dead = 0.02 N if needed.
ctrl.k_pos = k_pos;
ctrl.k_neg = k_neg;
ctrl.d_pos = d_pos;
ctrl.d_neg = d_neg;

%% 6) Save for Simulink
ctrl.A = A; ctrl.B = B; ctrl.C = C;
ctrl.K = K; ctrl.L = L;
% Encoder convention: theta_enc = 0 at upright (opposite zero from EOM).
% So phi (deviation from upright, used by the linearized model) is just
% theta_enc, wrapped to [-pi, pi]. Verify rotation direction matches the
% EOM before trusting K; flip sign of phi (and of omega_hat feedback) if
% it doesn't.
ctrl.theta_eq    = 0;          % subtract from encoder theta to get phi
ctrl.theta_sign  = +1;         % set to -1 if encoder rotation is flipped
ctrl.x_eq        = 0;
ctrl.Q_lqr = Q_lqr; ctrl.R_lqr = R_lqr;
ctrl.Q_kf  = Q_kf;  ctrl.R_kf  = R_kf;

out_file = fullfile('sysID','results', 'lqr_kalman.mat');
save(out_file, '-struct', 'ctrl');
fprintf('\nSaved controller + observer to %s\n', out_file);

%% 7) Quick sanity sim of the linear closed loop
t = 0:0.005:5;
z0 = [0; 0; 0.15; 0];          % 0.15 rad ~ 8.6 deg initial tilt
A_cl = [A - B*K,    B*K;
        zeros(4),   A - L*C];
% lifted state: [true state; estimation error]. Estimator gets correct y.
sys_cl = ss(A_cl, zeros(8,1), [C zeros(2,4)], 0);
[y_sim, ~, x_sim] = initial(sys_cl, [z0; z0], t);   % start with full err

if 0
figure('Name','Linear LQR+KF closed loop');
subplot(2,1,1); plot(t, y_sim(:,1)); ylabel('x [m]'); grid on;
title('Linear closed-loop response from initial tilt');
subplot(2,1,2); plot(t, y_sim(:,2)); ylabel('\phi [rad]'); xlabel('t [s]'); grid on;
end

%% ---- local function: EOM with force as input (theta = pi + phi) ----
function dz = eom_force(z, F, M, m, b, c, l)
    g = 9.81;
    x_v   = z(2);
    phi   = z(3);
    omega = z(4);
    theta = pi + phi;            % map back to original angle convention

    Delta = M + m*sin(theta)^2;
    num1  = F - c*x_v - m*l*omega^2*sin(theta) ...
            - m*g*cos(theta)*sin(theta) - (b/l)*cos(theta)*omega;
    num2  = m*l*cos(theta)*num1 - (M+m)*(m*g*l*sin(theta) + b*omega);

    dz = zeros(4,1);
    dz(1) = x_v;
    dz(2) = num1 / Delta;
    dz(3) = omega;
    dz(4) = num2 / (m*l^2 * Delta);
end
