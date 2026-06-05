function plot_results(x_mpc, u_mpc, x_lqr, u_lqr, param, mode)
fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

Ts = param.Ts;
T_mpc = size(u_mpc, 2);
T_lqr = size(u_lqr, 2);
t_mpc = (0:T_mpc) * Ts;
t_lqr = (0:T_lqr) * Ts;

%% Figure 1: 2D trajectory (x-y plane)
figure('Name', 'Trajectory', 'Position', [100 100 600 500]);
plot(x_mpc(2,:), x_mpc(1,:), 'b-', 'LineWidth', 1.5); hold on;
plot(x_lqr(2,:), x_lqr(1,:), 'r--', 'LineWidth', 1.5);
plot(x_mpc(2,1), x_mpc(1,1), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
plot(0, 0, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

xlabel('\delta y [m]'); ylabel('\delta x [m]');
title('Relative Trajectory (Hill Frame)');
grid minor; axis equal;

legend('MPC', 'LQR', 'Start', 'Target');
exportgraphics(gcf, fullfile(fig_dir, ['trajectory_' mode '.pdf']), 'ContentType', 'vector');

%% Figure 2: State trajectories over time
figure('Name', 'States', 'Position', [100 100 800 600]);

subplot(2,2,1)
plot(t_mpc, x_mpc(1,:), 'b-', 'LineWidth', 1.5); hold on;
plot(t_lqr, x_lqr(1,:), 'r--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('\delta x [m]');
title('Radial Position'); grid minor;
legend('MPC', 'LQR', 'Location', 'best');

subplot(2,2,2)
plot(t_mpc, x_mpc(2,:), 'b-', 'LineWidth', 1.5); hold on;
plot(t_lqr, x_lqr(2,:), 'r--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('\delta y [m]');
title('In-track Position'); grid minor;

subplot(2,2,3)
plot(t_mpc, x_mpc(3,:), 'b-', 'LineWidth', 1.5); hold on;
plot(t_lqr, x_lqr(3,:), 'r--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('d(dx)/dt [m/s]');
title('Radial Velocity'); grid minor;

subplot(2,2,4)
plot(t_mpc, x_mpc(4,:), 'b-', 'LineWidth', 1.5); hold on;
plot(t_lqr, x_lqr(4,:), 'r--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('d(dy)/dt [m/s]');
title('In-track Velocity'); grid minor;

exportgraphics(gcf, fullfile(fig_dir, ['states_' mode '.pdf']), 'ContentType', 'vector');

%% Figure 3: Control inputs over time
figure('Name', 'Inputs', 'Position', [100 100 800 400]);

subplot(1,2,1)
stairs((0:T_mpc-1)*Ts, u_mpc(1,:), 'b-', 'LineWidth', 1.5); hold on;
stairs((0:T_lqr-1)*Ts, u_lqr(1,:), 'r--', 'LineWidth', 1.5);
yline(param.u_max, 'k:', 'LineWidth', 1);
yline(-param.u_max, 'k:', 'LineWidth', 1);
xlabel('Time [s]'); ylabel('u_x [m/s^2]');
title('Radial Thrust'); grid minor;
legend('MPC', 'LQR', 'u_{max}', 'Location', 'best');

subplot(1,2,2)
stairs((0:T_mpc-1)*Ts, u_mpc(2,:), 'b-', 'LineWidth', 1.5); hold on;
stairs((0:T_lqr-1)*Ts, u_lqr(2,:), 'r--', 'LineWidth', 1.5);
yline(param.u_max, 'k:', 'LineWidth', 1);
yline(-param.u_max, 'k:', 'LineWidth', 1);
xlabel('Time [s]'); ylabel('u_y [m/s^2]');
title('In-track Thrust'); grid minor;

exportgraphics(gcf, fullfile(fig_dir, ['inputs_' mode '.pdf']), 'ContentType', 'vector');

%% Figure 4: Thrust magnitude
figure('Name', 'Thrust magnitude', 'Position', [100 100 600 300]);
u_mag_mpc = sqrt(u_mpc(1,:).^2 + u_mpc(2,:).^2);
u_mag_lqr = sqrt(u_lqr(1,:).^2 + u_lqr(2,:).^2);
stairs((0:T_mpc-1)*Ts, u_mag_mpc, 'b-', 'LineWidth', 1.5); hold on;
stairs((0:T_lqr-1)*Ts, u_mag_lqr, 'r--', 'LineWidth', 1.5);
yline(param.u_max, 'k:', 'LineWidth', 1);
xlabel('Time [s]'); ylabel('||u|| [m/s^2]');
title('Thrust Magnitude'); grid minor;
legend('MPC', 'LQR', 'u_{max}', 'Location', 'best');

exportgraphics(gcf, fullfile(fig_dir, ['thrust_mag_' mode '.pdf']), 'ContentType', 'vector');

end
