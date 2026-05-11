clear
clc

shall_plot = 1;
N = 500;
skip = 0;
% list of data files to use for ID
data_files = { ...
    'real_data\chirp1_20260508_200223.mat', ...
    'real_data\chirp2_20260508_200322.mat', ...
    'real_data\1sin1_20260506_191241.mat', ...
    'real_data\1sin2_20260506_190731.mat', ...
    % 'real_data\1sin4_20260506_191140.mat', ...
    % 'real_data\multisin_20260506_191743.mat', ...
    % 'real_data\1_05_Doub_20260506_190103.mat', ...
    % 'real_data\05_1_Doub_20260506_185819.mat', ...
};

% create id data objects + remember each file's measured initial state
n_data = numel(data_files);
data_set = cell(1, n_data);
x0_mat = zeros(4, n_data); % rows: states; cols: experiments
x_all = []; theta_all = []; % pooled signals for OutputWeight
for i = 1:n_data
    S = load(data_files{i});
    x_vec = [S.x(1+skip:N+skip), S.theta(1+skip:N+skip)];
    outu_i = S.outu(1+skip:N+skip);
    data_set{i} = iddata(x_vec, outu_i, (S.t(2) - S.t(1)));
    x0_mat(:,i) = [S.x(1+skip); 0; S.theta(1+skip); 0]; % measured pos, zero velocity guess
    x_all     = [x_all;     S.x(1+skip:N+skip)];     %#ok<AGROW>
    theta_all = [theta_all; S.theta(1+skip:N+skip)]; %#ok<AGROW>
end
data_merged = merge(data_set{:});

% creating the object for parameter estimation
order = [2 1 4]; % [output input states]
%         [M,   m,      b, c,    l,  k_pos, k_neg, d_pos, d_neg,  Fc,    Ft]
params0 = [1, 0.1, 0.0005, 1, 0.38,    3,     3,    0.05,  0.05, 0.5, -0.03]; % initial guess

% per-experiment initial state: pass Nx-by-Ne matrix straight to the constructor
sys0 = idnlgrey('sysEOM_asym', order, params0, x0_mat, 0);

sys0.Parameters(1).Minimum = 0.5;  sys0.Parameters(1).Maximum = 1;
sys0.Parameters(2).Minimum = 0;    sys0.Parameters(2).Maximum = 0.3;
sys0.Parameters(3).Minimum = 0;    sys0.Parameters(3).Maximum = 0.02;
sys0.Parameters(4).Minimum = 0;    sys0.Parameters(4).Maximum = 2;
sys0.Parameters(5).Minimum = 0.35; sys0.Parameters(5).Maximum = 0.45;
sys0.Parameters(6).Minimum = 1;    sys0.Parameters(6).Maximum = 5;
sys0.Parameters(7).Minimum = 1;    sys0.Parameters(7).Maximum = 5;
sys0.Parameters(8).Minimum = 0;    sys0.Parameters(8).Maximum = 0.3;
sys0.Parameters(9).Minimum = 0;    sys0.Parameters(9).Maximum = 0.3;
sys0.Parameters(10).Minimum = 0;   sys0.Parameters(10).Maximum = 2;
sys0.Parameters(11).Minimum = -0.3;  sys0.Parameters(11).Maximum = 0.3;

% positions fixed to measured values, velocities free
sys0.InitialStates(1).Fixed = true;
sys0.InitialStates(2).Fixed = false;
sys0.InitialStates(3).Fixed = true;
sys0.InitialStates(4).Fixed = false;

opt = nlgreyestOptions;
opt.SearchMethod = 'auto';
opt.SearchOptions.MaxIterations = 50;
opt.OutputWeight = diag([1/var(x_all), 1/var(theta_all)]);
opt.Display = 'on';

% actual optimization
sys_id = nlgreyest(data_merged, sys0, opt);

sys_id_vec = getpvec(sys_id);

% per-experiment initial state recovered from the multi-experiment sys_id
x0_post = zeros(4, n_data);
for i = 1:n_data
    for k = 1:4
        v = sys_id.InitialStates(k).Value;
        x0_post(k,i) = v(i);
    end
end

% result visualization and saving (build a single-experiment model per file)
fit_mat     = zeros(2, n_data); % rows: [x; theta] fit %, cols: experiments
y_model_set = cell(1, n_data);
err_set     = cell(1, n_data);
for i = 1:n_data
    sys_i = idnlgrey('sysEOM_asym', order, sys_id_vec, x0_post(:,i), 0);
    [y_model, fit, ~] = compare(data_set{i}, sys_i);
    fit_mat(:,i)   = fit(:);
    y_model_set{i} = y_model.OutputData;          % [x_model, theta_model]
    err_set{i}     = data_set{i}.OutputData - y_model_set{i}; % real - model
end
disp('the worst fit to the data is: ')
disp(min(fit_mat(:)))

param.M     = sys_id_vec(1);
param.m     = sys_id_vec(2);
param.b     = sys_id_vec(3);
param.c     = sys_id_vec(4);
param.l     = sys_id_vec(5);
param.k_pos = sys_id_vec(6);
param.k_neg = sys_id_vec(7);
param.d_pos = sys_id_vec(8);
param.d_neg = sys_id_vec(9);
param.Fc    = sys_id_vec(10);
param.Ft = sys_id_vec(11);
disp(param);

% per-file plots: input, x (real vs model), theta (real vs model), errors
if shall_plot
for i = 1:n_data
    t_i        = (0:size(data_set{i}.OutputData,1)-1)' * data_set{i}.Ts;
    u_i        = data_set{i}.InputData;
    y_real     = data_set{i}.OutputData;
    y_mod      = y_model_set{i};
    e_i        = err_set{i};

    figure('Name', data_files{i});
    subplot(4,1,1);
    plot(t_i, u_i, 'k'); grid on;
    ylabel('u'); title(data_files{i}, 'Interpreter', 'none');

    subplot(4,1,2);
    plot(t_i, y_real(:,1), 'b', t_i, y_mod(:,1), 'r--'); grid on;
    ylabel('x');
    legend('real', sprintf('model (fit: %.2f%%)', fit_mat(1,i)), 'Location','best');

    subplot(4,1,3);
    plot(t_i, y_real(:,2), 'b', t_i, y_mod(:,2), 'r--'); grid on;
    ylabel('\theta');
    legend('real', sprintf('model (fit: %.2f%%)', fit_mat(2,i)), 'Location','best');

    subplot(4,1,4);
    e_x_n     = e_i(:,1) / sqrt(var(x_all));
    e_theta_n = e_i(:,2) / sqrt(var(theta_all));
    plot(t_i, e_x_n, 'b', t_i, e_theta_n, 'r'); grid on;
    ylabel('error / std'); xlabel('t [s]');
    legend('x error','\theta error','Location','best');
end
end