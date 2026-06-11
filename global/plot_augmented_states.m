function plot_augmented_states(x_hat, dx_hat, theta_hat, dtheta_hat, h, T, results_dir, timestamp)
% PLOT_AUGMENTED_STATES  Plot all 5 LQI augmented-state trajectories.
%   States: x, v, theta, omega, xi (integral of position error).
%   Saves augmented_states[_<timestamp>].png to results_dir.
if nargin < 8 || isempty(timestamp)
    timestamp = '';
end
suffix = '';
if ~isempty(timestamp), suffix = ['_' timestamp]; end

t_vec = 0:h:T;
xi_est = cumsum(-x_hat) * h;   % xi_dot = -x_hat  (x_ref = 0)

fig_states = figure('Name', 'Augmented state trajectories');
% subplot(5,1,1); plot(t_vec, x_hat, t_vec,ref_signal);      ylabel('x [m]');           grid on;
subplot(5,1,1); plot(t_vec, x_hat);      ylabel('x [m]');           grid on;
title('Augmented state trajectories');
subplot(5,1,2); plot(t_vec, dx_hat);     ylabel('v [m/s]');         grid on;
subplot(5,1,3); plot(t_vec, theta_hat);  ylabel('\theta [rad]');    grid on;
subplot(5,1,4); plot(t_vec, dtheta_hat); ylabel('\omega [rad/s]');  grid on;
subplot(5,1,5); plot(t_vec, xi_est);     ylabel('\xi [m{\cdot}s]'); xlabel('t [s]'); grid on;

saveas(fig_states, fullfile(results_dir, ['augmented_states' suffix '.png']));
end
