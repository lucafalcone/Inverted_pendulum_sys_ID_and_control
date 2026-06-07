function param = define_params()
%% Orbital parameters
param.mu    = 3.986e14;         % Earth gravitational parameter [m^3/s^2]
param.Re    = 6.371e6;          % Earth radius [m]
param.h_orb = 500e3;            % Orbital altitude [m]
param.R0    = param.Re + param.h_orb;  % Orbital radius [m]
param.n     = sqrt(param.mu / param.R0^3); % Mean motion [rad/s]
param.T_orb = 2*pi / param.n;  % Orbital period [s]

%% Spacecraft parameters
param.m_chaser = 500;           % Chaser mass [kg]

%% State and input dimensions
param.nx = 4;                   % Number of states: [dx, dy, dx_dot, dy_dot]
param.nu = 2;                   % Number of inputs: [u_x, u_y]

%% Discretisation
param.Ts = 1;                   % Sampling period [s]

%% Input constraints
param.u_max = 0.1;              % Maximum thrust acceleration [m/s^2]

%% Line-of-Sight (LOS) cone constraints
param.gamma_los = 10 * pi/180;  % LOS cone half-angle [rad]
param.r_p       = 2.5;          % Target platform radius [m]
param.r_tol     = 0.5;          % Constraint relaxation tolerance [m]

%% Soft-docking constraints
param.eta  = 1;                 % Soft-docking shape parameter
param.beta = 0.25;              % Soft-docking offset parameter

%% MPC horizons
param.N  = 30;                  % Prediction horizon (steps)

%% Cost function weights
param.Q = diag([10, 10, 100, 10]);  % MPc weight matrix
param.R = eye(param.nu);         % Input weights

% to keep lqr under constraint x velocity must be low
param.Q_lqr = diag([10, 10, 8000, 100]);  % LQR state weight matrix (independent of MPC)
param.R_lqr = 1e5 * eye(param.nu);     % LQR input weight matrix

%% Simulation settings
param.T_sim   = 200;            % Simulation time [s]
param.x0      = [-50; -5; 0.5; -0.5]; % Initial state
param.x_dock  = 0.05;            % Docking distance threshold [m]
param.r_crash = 0.1;              % Crash detection proximity radius [m]
param.v_crash = 0.1;            % Maximum safe approach speed within crash zone [m/s]

%% Disturbance / noise parameters
param.w_std    = 1e-2;          % Process noise
param.v_std    = 1e-2;          % Measurement noise

end
