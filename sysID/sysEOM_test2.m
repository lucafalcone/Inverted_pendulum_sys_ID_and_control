% Parameters
param.M   = 1.0;
param.m   = 0.2;
param.g   = 9.81;
param.b   = 0.01;
param.c   = 0.01;
param.l   = 0.3;
param.k_m = 1.0;

dt_imp = 0.01;
F_imp  = 100;
u = @(t) F_imp * (t <= dt_imp);

x0    = [0; 0; pi; 0];   % hanging down
tspan = [0 10];
[t, x] = ode45(@(t,x) sysEOM(param, x, u(t)), tspan, x0);

% --- Animation ---
figure('Color','w');
ax = axes;
axis equal; grid on;
xlabel('x [m]'); ylabel('y [m]');
title('Cart-Pendulum');

l      = param.l;
cw     = 0.15;   % cart width
ch     = 0.08;   % cart height

for k = 1:3:length(t)   % skip frames for speed

    cart_x = x(k,1);
    theta  = x(k,3);   % 0 = upright, pi = hanging

    % Pendulum bob position
    bob_x = cart_x + l*sin(theta);
    bob_y =          l*cos(theta);   % y=0 is pivot height

    cla(ax);

    % Rail
    line([-2 2], [0 0], 'Color',[0.7 0.7 0.7], 'LineWidth', 1);

    % Cart (rectangle)
    rectangle('Position', [cart_x - cw/2, -ch/2, cw, ch], ...
        'FaceColor', [0.2 0.5 0.8], 'EdgeColor','k', 'LineWidth',1.5);

    % Pendulum rod
    line([cart_x, bob_x], [0, bob_y], 'Color','k', 'LineWidth', 2);

    % Bob
    viscircles([bob_x, bob_y], 0.03, 'Color', [0.85 0.2 0.2], 'LineWidth', 2);

    % Axes limits that follow the cart
    xlim([cart_x - 1.5,  cart_x + 1.5]);
    ylim([-l - 0.2,       l + 0.2]);

    title(sprintf('t = %.2f s', t(k)));
    drawnow;
end