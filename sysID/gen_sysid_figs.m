% Generate sysID training-fit and validation-fit figures for the report.
% Uses RK4 integration with segmented restarts (no System ID Toolbox needed).

clear; clc;

proj = '/home/luca/Documents/GitHub/Inverted_pendulum_sys_ID_and_control';
addpath(fullfile(proj, 'sysID', 'EOMs'));
img_dir = fullfile(proj, 'report', 'SC_Integration_Project', 'images');

S_p = load(fullfile(proj, 'global', 'param_64_175.mat'));
p = S_p.param;
Ts = 0.01;

%% Training data fit  (chirp2, June 8, first file)
S = load(fullfile(proj, 'sysID', 'real_data', 'chirp2_20260608_113815.mat'));
N = 500;
t_tr   = S.t(1:N);
u_tr   = S.outu(1:N);
x_tr   = S.x(1:N);
th_tr  = S.theta(1:N);

x_sim_tr = rk4_restart(x_tr, th_tr, u_tr, Ts, p, 20);

% NRMSE
fit_x  = nrmse(x_tr,  x_sim_tr(1,:)');
fit_th = nrmse(th_tr, x_sim_tr(3,:)');
fprintf('Training fit:  x = %.1f%%,  theta = %.1f%%\n', fit_x, fit_th);

fig1 = figure('Visible','off','Position',[100 100 800 520]);
subplot(3,1,1);
plot(t_tr, u_tr, 'k', 'LineWidth', 0.8); grid on;
ylabel('$u$ [V]','Interpreter','latex'); xlim([t_tr(1) t_tr(end)]);
title('Training data fit (chirp2, open-loop experiment)');

subplot(3,1,2);
plot(t_tr, x_tr, 'b', 'LineWidth', 1); hold on;
plot(t_tr, x_sim_tr(1,:)', 'r--', 'LineWidth', 1); grid on;
ylabel('$x$ [m]','Interpreter','latex'); xlim([t_tr(1) t_tr(end)]);
legend(sprintf('measured'), sprintf('model (%.0f%%)', fit_x), 'Location','best');

subplot(3,1,3);
plot(t_tr, th_tr, 'b', 'LineWidth', 1); hold on;
plot(t_tr, x_sim_tr(3,:)', 'r--', 'LineWidth', 1); grid on;
ylabel('$\theta$ [rad]','Interpreter','latex');
xlabel('$t$ [s]','Interpreter','latex'); xlim([t_tr(1) t_tr(end)]);
legend(sprintf('measured'), sprintf('model (%.0f%%)', fit_th), 'Location','best');

exportgraphics(fig1, fullfile(img_dir, 'sysid_training.png'), 'Resolution', 150);
fprintf('Saved sysid_training.png\n');

%% Validation data fit  (1sin1, skip first 500 samples)
S2 = load(fullfile(proj, 'sysID', 'real_data', '1sin1_20260506_191241.mat'));
skip = 500;
t_va   = S2.t(1+skip : N+skip);
u_va   = S2.outu(1+skip : N+skip);
x_va   = S2.x(1+skip : N+skip);
th_va  = S2.theta(1+skip : N+skip);

x_sim_va = rk4_restart(x_va, th_va, u_va, Ts, p, 20);

fit_x_v  = nrmse(x_va,  x_sim_va(1,:)');
fit_th_v = nrmse(th_va, x_sim_va(3,:)');
fprintf('Validation fit: x = %.1f%%,  theta = %.1f%%\n', fit_x_v, fit_th_v);

fig2 = figure('Visible','off','Position',[100 100 800 520]);
subplot(3,1,1);
plot(t_va, u_va, 'k', 'LineWidth', 0.8); grid on;
ylabel('$u$ [V]','Interpreter','latex'); xlim([t_va(1) t_va(end)]);
title('Validation data fit (sinusoidal, unseen data)');

subplot(3,1,2);
plot(t_va, x_va, 'b', 'LineWidth', 1); hold on;
plot(t_va, x_sim_va(1,:)', 'r--', 'LineWidth', 1); grid on;
ylabel('$x$ [m]','Interpreter','latex'); xlim([t_va(1) t_va(end)]);
legend(sprintf('measured'), sprintf('model (%.0f%%)', fit_x_v), 'Location','best');

subplot(3,1,3);
plot(t_va, th_va, 'b', 'LineWidth', 1); hold on;
plot(t_va, x_sim_va(3,:)', 'r--', 'LineWidth', 1); grid on;
ylabel('$\theta$ [rad]','Interpreter','latex');
xlabel('$t$ [s]','Interpreter','latex'); xlim([t_va(1) t_va(end)]);
legend(sprintf('measured'), sprintf('model (%.0f%%)', fit_th_v), 'Location','best');

exportgraphics(fig2, fullfile(img_dir, 'sysid_validation_new.png'), 'Resolution', 150);
fprintf('Saved sysid_validation_new.png\n');

% -------------------------------------------------------------------------
function x_sim = rk4_restart(x_meas, th_meas, u_vec, Ts, p, seg)
    N = length(u_vec);
    x_sim = zeros(4, N);
    for k = 1:N-1
        if k == 1 || mod(k-1, seg) == 0
            v0 = 0;  w0 = 0;
            if k > 1
                v0 = (x_meas(k) - x_meas(k-1)) / Ts;
                w0 = (th_meas(k) - th_meas(k-1)) / Ts;
            end
            x_sim(:,k) = [x_meas(k); v0; th_meas(k); w0];
        end
        xk = x_sim(:,k);
        uk = u_vec(k);
        k1 = f_eom(xk,       uk, Ts, p);
        k2 = f_eom(xk+Ts/2*k1, uk, Ts, p);
        k3 = f_eom(xk+Ts/2*k2, uk, Ts, p);
        k4 = f_eom(xk+Ts*k3,   uk, Ts, p);
        x_sim(:,k+1) = xk + Ts/6*(k1+2*k2+2*k3+k4);
    end
    x_sim(:,N) = [x_meas(N); 0; th_meas(N); 0];
end

function dx = f_eom(x, u, ~, p)
    [dx_raw, ~] = sysEOM_asym(0, x, u, p.M, p.m, p.b, p.c, p.l, ...
                               p.k_pos, p.k_neg, p.d_pos, p.d_neg, p.Fc);
    dx = dx_raw(:);
end

function fit = nrmse(y_real, y_model)
    err = y_real - y_model;
    fit = max(0, 1 - norm(err)/norm(y_real - mean(y_real))) * 100;
end
