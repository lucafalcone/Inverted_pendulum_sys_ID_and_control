function plot_kalman_comparison(x, theta, x_hat, theta_hat, dx_hat, dtheta_hat, h, T, results_dir, timestamp)
% PLOT_KALMAN_COMPARISON  Plot measured vs Kalman-estimated states.
%   Saves kalman_comparison[_<timestamp>].png to results_dir.
if nargin < 10 || isempty(timestamp)
    timestamp = '';
end
suffix = '';
if ~isempty(timestamp), suffix = ['_' timestamp]; end

t_vec = 0:h:T;

fig_kf = figure('Name', 'Kalman filter comparison');

subplot(4,1,1);
plot(t_vec, x, t_vec, x_hat); ylabel('x [m]'); grid on;
legend('measured', 'Kalman'); title('Kalman filter comparison');

subplot(4,1,2);
plot(t_vec, theta, t_vec, theta_hat); ylabel('\theta [rad]'); grid on;

subplot(4,1,3);
plot(0:h:T-h, diff(x)./h, t_vec, dx_hat); ylabel('v [m/s]'); grid on;

subplot(4,1,4);
plot(0:h:T-h, diff(theta)./h, t_vec, dtheta_hat);
ylabel('\omega [rad/s]'); xlabel('t [s]'); grid on;

saveas(fig_kf, fullfile(results_dir, ['kalman_comparison' suffix '.png']));
end
