clear
clc

N = 250;
skip = 0;
% list of data files to use for ID
data_files = { ...
    'real_data\chirp1_20260508_200223.mat', ...
    'real_data\chirp2_20260508_200322.mat', ...
    'real_data\1sin1_20260506_191241.mat', ...
    'real_data\1sin2_20260506_190731.mat', ...
    'real_data\1sin4_20260506_191140.mat', ...
    'real_data\multisin_20260506_191743.mat', ...
};

% create id data objects + remember each file's measured initial state
n_data = numel(data_files);
data_set = cell(1, n_data);
x0_mat = zeros(4, n_data); % rows: states; cols: experiments
for i = 1:n_data
    S = load(data_files{i});
    x_vec = [S.x(1+skip:N+skip), S.theta(1+skip:N+skip)];
    outu_i = S.outu(1+skip:N+skip);
    data_set{i} = iddata(x_vec, outu_i, (S.t(2) - S.t(1)));
    x0_mat(:,i) = [S.x(1+skip); 0; S.theta(1+skip); 0]; % measured pos, zero velocity guess
end
% keep last loaded x/theta around for the OutputWeight below
x = S.x; theta = S.theta;
data_merged = merge(data_set{:});

% creating the object for parameter estimation
order = [2 1 4]; % [output input states]
%         [  M,   m,      b, c,    l,   k]
params0 = [0.8, 0.1, 0.0005, 1, 0.38,   3]; % initial guess

% per-experiment initial state: pass Nx-by-Ne matrix straight to the constructor
sys0 = idnlgrey('sysEOM', order, params0, x0_mat, 0);

sys0.Parameters(1).Minimum = 0;    % sys0.Parameters(1).Maximum = 3;
sys0.Parameters(2).Minimum = 0;    % sys0.Parameters(2).Maximum = 1;
sys0.Parameters(3).Minimum = 0;    % sys0.Parameters(3).Maximum = 0.5;
sys0.Parameters(4).Minimum = 0;    % sys0.Parameters(4).Maximum = 10;
sys0.Parameters(5).Minimum = 0.3;  % sys0.Parameters(5).Maximum = 0.45;
sys0.Parameters(6).Minimum = 0;    % sys0.Parameters(6).Maximum = 8;

% positions fixed to measured values, velocities free
sys0.InitialStates(1).Fixed = true;
sys0.InitialStates(2).Fixed = false;
sys0.InitialStates(3).Fixed = true;
sys0.InitialStates(4).Fixed = false;

opt = nlgreyestOptions;
opt.SearchMethod = 'auto';
opt.SearchOptions.MaxIterations = 50;
opt.OutputWeight = diag([1/var(x), 1/var(theta)]);
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
fit_vec = [];
for i = 1:n_data
    sys_i = idnlgrey('sysEOM', order, sys_id_vec, x0_post(:,i), 0);
    [~, fit, ~] = compare(data_set{i}, sys_i);
    fit_vec = [fit_vec; fit]; %#ok<AGROW>
end
disp('the worst fit to the data is: ')
disp(min(fit_vec))

param.M = sys_id_vec(1);
param.m = sys_id_vec(2);
param.b = sys_id_vec(3);
param.c = sys_id_vec(4);
param.l = sys_id_vec(5);
param.k = sys_id_vec(6);
disp(param);

for i = 1:n_data
    sys_i = idnlgrey('sysEOM', order, sys_id_vec, x0_post(:,i), 0);
    figure;
    compare(data_set{i}, sys_i);
    title(data_files{i}, 'Interpreter', 'none');
end

%% SAVE OUTPUT
save('param.mat','param');