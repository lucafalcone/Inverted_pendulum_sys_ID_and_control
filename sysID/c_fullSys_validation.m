clear
clc

N = 500;
% list of data files to use for validation (should NOT overlap with ID set)
data_files = { 
    'real_data\1sin1_20260506_191241.mat', ...
    'real_data\1sin2_20260506_190731.mat', ...
    'real_data\1sin4_20260506_191140.mat', ...
    'real_data\multisin_20260506_191743.mat', ...
};

% load identified parameters (expects struct `param` with M, m, b, c, l, k)
load('param.mat');
params_id = [param.M, param.m, param.b, param.c, param.l, param.k];

% create id data objects + remember each file's measured initial state
n_data = numel(data_files);
data_set = cell(1, n_data);
x0_set = cell(1, n_data);
for i = 1:n_data
    S = load(data_files{i});
    x_vec = [S.x(1:N), S.theta(1:N)];
    outu_i = S.outu(1:N);
    data_set{i} = iddata(x_vec, outu_i, (S.t(2) - S.t(1)));
    x0_set{i} = [S.x(1); 0; S.theta(1); 0]; % measured pos, zero velocity guess
end

% build the model with the identified parameters (no estimation)
order = [2 1 4]; % [output input states]
sys_id = idnlgrey('sysEOM', order, params_id, zeros(4,1), 0); % placeholder x0; overwritten per experiment below

% strict open-loop validation: lock initial states to measured values
sys_id.InitialStates(1).Fixed = true;
sys_id.InitialStates(2).Fixed = true;
sys_id.InitialStates(3).Fixed = true;
sys_id.InitialStates(4).Fixed = true;

% fits on each validation set (re-set initial state per file)
fit_vec = [];
for i = 1:n_data
    for k = 1:4
        sys_id.InitialStates(k).Value = x0_set{i}(k);
    end
    [~, fit, ~] = compare(data_set{i}, sys_id);
    fit_vec = [fit_vec; fit]; %#ok<AGROW>
end
disp('the average fit on validation data is: ')
disp(mean(fit_vec))

for i = 1:n_data
    for k = 1:4
        sys_id.InitialStates(k).Value = x0_set{i}(k);
    end
    figure;
    compare(data_set{i}, sys_id);
    title(data_files{i}, 'Interpreter', 'none');
end