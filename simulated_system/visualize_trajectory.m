function visualize_trajectory(t, X, param, opts)
% VISUALIZE_TRAJECTORY  Animate cart-pendulum and trace both trajectories.
%
%   visualize_trajectory(t, X, param)
%   visualize_trajectory(t, X, param, opts)
%
%   t      : Nx1 time vector
%   X      : Nx4 state matrix, columns = [s, s_dot, theta, theta_dot]
%            theta measured from upward vertical (theta=0 -> pendulum up)
%   param  : struct with field l (pendulum length)
%   opts   : optional struct
%            .speed     playback speed multiplier (default 1)
%            .save_gif  filename to write an animated gif (default '')
%            .trail     seconds of trail to keep (default inf = full)

    arguments
        t      (:,1) double
        X      (:,4) double
        param  struct
        opts.speed     (1,1) double = 1
        opts.save_gif  (1,:) char   = ''
        opts.trail     (1,1) double = inf
    end

    l = param.l;

    s     = X(:,1);
    theta = X(:,3);

    % bob position in world frame (theta = 0 -> straight up)
    bob_x = s + l*sin(theta);
    bob_y =     l*cos(theta);

    % figure styling
    fig = figure('Color','w','Position',[200 200 900 520]);
    ax  = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    ax.GridAlpha   = 0.15;
    ax.FontName    = 'Helvetica';
    ax.FontSize    = 11;
    ax.LineWidth   = 1;
    ax.DataAspectRatio = [1 1 1];

    pad   = 0.3;
    xlim(ax, [min([s; bob_x])-pad, max([s; bob_x])+pad]);
    ylim(ax, [-1.3*l, 1.5*l]);
    xlabel(ax,'x [m]'); ylabel(ax,'y [m]');
    title(ax,'Cart-pendulum trajectory','FontWeight','normal');

    % ground / track
    yline(ax, 0, 'Color',[0.3 0.3 0.3], 'LineWidth', 1.2);

    % colors
    cartColor = [0.20 0.45 0.85];
    rodColor  = [0.15 0.15 0.20];
    bobColor  = [0.85 0.30 0.25];
    trailCart = [0.20 0.45 0.85 0.35];
    trailBob  = [0.85 0.30 0.25 0.35];

    % cart geometry
    cartW = 0.25; cartH = 0.14;

    % graphics handles (initialised at first frame)
    hTrailCart = plot(ax, nan, nan, '-', 'Color', trailCart, 'LineWidth', 1.5);
    hTrailBob  = plot(ax, nan, nan, '-', 'Color', trailBob,  'LineWidth', 1.5);
    hCart      = rectangle(ax, 'Position',[s(1)-cartW/2, -cartH/2, cartW, cartH], ...
                           'Curvature',0.15, 'FaceColor', cartColor, ...
                           'EdgeColor','k', 'LineWidth',1);
    hRod       = plot(ax, [s(1), bob_x(1)], [0, bob_y(1)], '-', ...
                      'Color', rodColor, 'LineWidth', 2.5);
    hBob       = plot(ax, bob_x(1), bob_y(1), 'o', ...
                      'MarkerSize', 11, 'MarkerFaceColor', bobColor, ...
                      'MarkerEdgeColor','k', 'LineWidth',1);
    hPivot     = plot(ax, s(1), 0, 'o', ...
                      'MarkerSize', 5, 'MarkerFaceColor','k', ...
                      'MarkerEdgeColor','k');
    hTime      = text(ax, 0.02, 0.96, '', 'Units','normalized', ...
                      'FontName','Helvetica', 'FontSize', 11, ...
                      'VerticalAlignment','top', ...
                      'BackgroundColor',[1 1 1 0.7]);

    legend(ax, [hTrailCart, hTrailBob], {'cart path','bob path'}, ...
           'Location','northeast','Box','off');

    N = numel(t);

    if ~isempty(opts.save_gif)
        % gif export: write every sample, pace via DelayTime
        dt = [diff(t); mean(diff(t))];
        for k = 1:N
            updateFrame(k);
            drawnow;
            frame    = getframe(fig);
            [A, map] = rgb2ind(frame2im(frame), 256);
            if k == 1
                imwrite(A, map, opts.save_gif, 'gif', ...
                        'LoopCount', inf, 'DelayTime', dt(k)/opts.speed);
            else
                imwrite(A, map, opts.save_gif, 'gif', ...
                        'WriteMode','append', 'DelayTime', dt(k)/opts.speed);
            end
        end
        return
    end

    % real-time playback: index follows wall-clock, frames drop if slow
    targetFPS = 60;
    frameDt   = 1/targetFPS;
    t0        = t(1);
    tStart    = tic;
    nextFrame = frameDt;
    k         = 1;
    while ishandle(fig) && k < N
        elapsed = toc(tStart) * opts.speed;
        simT    = t0 + elapsed;
        while k < N && t(k+1) <= simT
            k = k + 1;
        end
        updateFrame(k);
        drawnow limitrate;
        wait = (nextFrame - toc(tStart)) / opts.speed;
        if wait > 0, pause(wait); end
        nextFrame = nextFrame + frameDt;
    end
    if ishandle(fig), updateFrame(N); drawnow; end

    function updateFrame(k)
        if isfinite(opts.trail)
            k0 = find(t >= t(k) - opts.trail, 1, 'first');
        else
            k0 = 1;
        end
        set(hTrailCart, 'XData', s(k0:k),     'YData', zeros(k-k0+1,1));
        set(hTrailBob,  'XData', bob_x(k0:k), 'YData', bob_y(k0:k));
        set(hCart, 'Position', [s(k)-cartW/2, -cartH/2, cartW, cartH]);
        set(hRod,  'XData', [s(k), bob_x(k)], 'YData', [0, bob_y(k)]);
        set(hBob,  'XData', bob_x(k), 'YData', bob_y(k));
        set(hPivot,'XData', s(k));
        set(hTime, 'String', sprintf(' t = %5.2f s ', t(k)));
    end
end
