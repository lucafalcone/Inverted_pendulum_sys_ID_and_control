%%
clear all
clc
close all

set(0,'DefaultFigureVisible','on') 

script_dir = fileparts(mfilename('fullpath'));

%% 1. Load parameters and build model
param = define_params();
LTI   = build_cw_model(param);

%% 2. Set up dimensions
dim.nx = param.nx;
dim.nu = param.nu;
dim.N  = param.N;

%% 3. Compute terminal ingredients
fprintf('=== Computing terminal ingredients ===\n');
terminal = compute_terminal(LTI, param);

%% 4. Build prediction model and cost function
weight.Q = param.Q;
weight.R = param.R;
weight.P = terminal.P;

predmod = predmodgen(LTI, dim);
[H, h, ~] = costgen(predmod, weight, dim);

%% 5. Build constraint matrices
constraints = constraintgen(LTI, param, dim);

%% 6. Run MPC closed-loop simulation
T_sim = param.T_sim;
nx    = dim.nx;
nu    = dim.nu;

x_mpc = zeros(nx, T_sim+1);
u_mpc = zeros(nu, T_sim);
J_mpc = zeros(1, T_sim);       % Cost at each step
t_solve = zeros(1, T_sim);     % Computation time per step

x_mpc(:,1) = param.x0;

fprintf('=== Running MPC simulation ===\n');

crashed = false;
for k = 1:T_sim
    x_0 = x_mpc(:,k);

    % Solve QP
    tic;
    [u_opt, fval, exitflag] = mpc_solve(x_0, H, h, constraints, dim);
    t_solve(k) = toc;

    % Apply first input only (receding horizon)
    u_mpc(:,k) = u_opt(1:nu);

    % Store cost
    J_mpc(k) = fval;

    % Propagate state
    x_mpc(:,k+1) = LTI.A * x_0 + LTI.B * u_mpc(:,k);

    pos_norm = norm(x_mpc(1:2, k+1));
    vel_norm = norm(x_mpc(3:4, k+1));

    % Check crash condition (within proximity zone at unsafe speed)
    if pos_norm < param.r_crash && vel_norm > param.v_crash
        fprintf('*** CRASH DETECTED at step k = %d (t = %.1f s): pos = %.2f m, vel = %.4f m/s ***\n', ...
            k, k*param.Ts, pos_norm, vel_norm);
        x_mpc = x_mpc(:, 1:k+1);
        u_mpc = u_mpc(:, 1:k);
        J_mpc = J_mpc(1:k);
        t_solve = t_solve(1:k);
        crashed = true;
        break;
    end

    % Check docking condition
    if pos_norm < param.x_dock
        fprintf('Docking achieved at step k = %d (t = %.1f s)\n', k, k*param.Ts);
        % Truncate trajectories
        x_mpc = x_mpc(:, 1:k+1);
        u_mpc = u_mpc(:, 1:k);
        J_mpc = J_mpc(1:k);
        t_solve = t_solve(1:k);
        break;
    end
end

T_actual = size(u_mpc, 2);
fprintf('Final position: [%.4f, %.4f] m\n', x_mpc(1,end), x_mpc(2,end));
fprintf('Final velocity: [%.4f, %.4f] m/s\n', x_mpc(3,end), x_mpc(4,end));
fprintf('Mean solve time: %.4f ms\n', mean(t_solve)*1000);

%% 7. Run LQR baseline for comparison
fprintf('\n=== Running LQR baseline ===\n');
[x_lqr, u_lqr, u_lqr_cmd] = lqr_baseline(LTI, param);

%% 8. Compute performance metrics
% Total delta-v (fuel consumption)
dv_mpc = sum(sqrt(u_mpc(1,:).^2 + u_mpc(2,:).^2)) * param.Ts;
dv_lqr = sum(sqrt(u_lqr(1,:).^2 + u_lqr(2,:).^2)) * param.Ts;

fprintf('\n=== Performance Comparison ===\n');
fprintf('MPC total delta-v:  %.4f m/s\n', dv_mpc);
fprintf('LQR total delta-v:  %.4f m/s\n', dv_lqr);

%% 9. Generate plots
plot_results(x_mpc, u_mpc, x_lqr, u_lqr, param, 'baseline');

%% 10. Save results
results.x_mpc    = x_mpc;
results.u_mpc    = u_mpc;
results.x_lqr     = x_lqr;
results.u_lqr     = u_lqr;
results.u_lqr_cmd = u_lqr_cmd;
results.J_mpc    = J_mpc;
results.t_solve  = t_solve;
results.dv_mpc   = dv_mpc;
results.dv_lqr   = dv_lqr;
results.param    = param;
results.terminal = terminal;

save(fullfile(script_dir, 'results_baseline.mat'), 'results');
fprintf('\nResults saved to results_baseline.mat\n');
