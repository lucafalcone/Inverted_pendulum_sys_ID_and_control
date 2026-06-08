function plot_kalman_comparison(x, theta, x_hat, theta_hat, dx_hat, dtheta_hat, h, T, results_dir)
% PLOT_KALMAN_COMPARISON  Plot and save Kalman filter comparison figures.
%   Produces two figures:
%     kalman_comparison.png  – measured vs estimated x, theta, v, omega
%     augmented_states.png   – all 5 augmented-state trajectories

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

saveas(fig_kf, fullfile(results_dir, 'kalman_comparison.png'));

% All 5 augmented states: x, v, theta, omega, xi
xi_est = cumsum(-x_hat) * h;   % xi_dot = x_ref - x_hat = -x_hat (x_ref = 0)

fig_states = figure('Name', 'Augmented state trajectories');
subplot(5,1,1); plot(t_vec, x_hat);      ylabel('x [m]');             grid on;
title('Augmented state trajectories');
subplot(5,1,2); plot(t_vec, dx_hat);     ylabel('v [m/s]');           grid on;
subplot(5,1,3); plot(t_vec, theta_hat);  ylabel('\theta [rad]');      grid on;
subplot(5,1,4); plot(t_vec, dtheta_hat); ylabel('\omega [rad/s]');    grid on;
subplot(5,1,5); plot(t_vec, xi_est);     ylabel('\xi [m{\cdot}s]'); xlabel('t [s]'); grid on;

saveas(fig_states, fullfile(results_dir, 'augmented_states.png'));
end
