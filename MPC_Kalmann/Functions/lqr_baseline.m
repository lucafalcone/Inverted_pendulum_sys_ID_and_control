function [x_lqr, u_lqr, u_lqr_cmd] = lqr_baseline(LTI, param)
T_sim = param.T_sim;
nx    = param.nx;
nu    = param.nu;

%% Compute LQR gain (uses R_lqr, tuned separately from MPC's R)
[~, ~, Kd] = dare(LTI.A, LTI.B, param.Q_lqr, param.R_lqr);
K = -Kd;

%% Simulate closed-loop system
x_lqr     = zeros(nx, T_sim+1);
u_lqr     = zeros(nu, T_sim);
u_lqr_cmd = zeros(nu, T_sim);
x_lqr(:,1) = param.x0;

for k = 1:T_sim
    % Compute LQR control input (unconstrained)
    u_k = K * x_lqr(:,k);
    u_lqr_cmd(:,k) = u_k;

    % Apply input saturation (post-hoc clipping)
    u_k = max(-param.u_max, min(param.u_max, u_k));

    u_lqr(:,k) = u_k;

    % Propagate state
    x_lqr(:,k+1) = LTI.A * x_lqr(:,k) + LTI.B * u_k;

    pos_norm = norm(x_lqr(1:2, k+1));
    vel_norm = norm(x_lqr(3:4, k+1));

    if pos_norm < param.r_crash
        disp("nearby sat");
    end

    % Crash detection: inside proximity zone at unsafe speed
    if pos_norm < param.r_crash && vel_norm > param.v_crash
        fprintf('*** LQR CRASH DETECTED at step k = %d (t = %.1f s): pos = %.2f m, vel = %.4f m/s ***\n', ...
            k, k*param.Ts, pos_norm, vel_norm);
        x_lqr     = x_lqr(:, 1:k+1);
        u_lqr     = u_lqr(:, 1:k);
        u_lqr_cmd = u_lqr_cmd(:, 1:k);
        return;
    end

    % Docking condition: within threshold at safe speed
    if pos_norm < param.x_dock && vel_norm <= param.v_crash
        fprintf('LQR docking achieved at step k = %d (t = %.1f s)\n', k, k*param.Ts);
        x_lqr     = x_lqr(:, 1:k+1);
        u_lqr     = u_lqr(:, 1:k);
        u_lqr_cmd = u_lqr_cmd(:, 1:k);
        return;
    end
end

end
