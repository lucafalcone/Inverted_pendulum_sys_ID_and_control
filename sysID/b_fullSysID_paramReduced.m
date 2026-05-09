clear
clc

N = 300;
% create id data object using real data 
load real_data\1_05_Doub_20260506_190103.mat
x_vec = [x(1:N), theta(1:N)]; outu = outu(1:N);
data1 = iddata(x_vec, outu, (t(2) - t(1)));
load real_data\05_1_Doub_20260506_185819.mat
x_vec = [x(1:N), theta(1:N)]; outu = outu(1:N);
data2 = iddata(x_vec, outu, (t(2) - t(1)));
load real_data\1sin1_20260506_191241.mat
x_vec = [x(1:N), theta(1:N)]; outu = outu(1:N);
data3 = iddata(x_vec, outu, (t(2) - t(1)));
load real_data\1sin2_20260506_190731.mat
x_vec = [x(1:N), theta(1:N)]; outu = outu(1:N);
data4 = iddata(x_vec, outu, (t(2) - t(1)));
data_merged = merge(data1, data2, data3, data4);

% creating the object for parameter estimation
order = [2 1 4]; % [output input states]
%         [  M,   m,   c, k_m]
params0 = [0.4, 0.1, 0.1,   3]; % initial guess
x0 = [x(1) 0 theta(1) 0]'; % fixed or initial guess if initial state not fixed

sys0 = idnlgrey('sysEOM_paramReduced', order, params0, x0, 0);

sys0.Parameters(1).Minimum = 0;     sys0.Parameters(1).Maximum = 1;
sys0.Parameters(2).Minimum = 0;     sys0.Parameters(2).Maximum = 1;
sys0.Parameters(3).Minimum = 0;     sys0.Parameters(3).Maximum = 0.5;
sys0.Parameters(4).Minimum = 0;     sys0.Parameters(4).Maximum = 10;
sys0.InitialStates(1).Fixed = false;
sys0.InitialStates(2).Fixed = false;
sys0.InitialStates(3).Fixed = false;
sys0.InitialStates(4).Fixed = false;

opt = nlgreyestOptions;
opt.SearchMethod = 'auto';
opt.SearchOptions.MaxIterations = 50;
opt.OutputWeight = diag([1/var(x), 1/var(theta)]);
opt.Display = 'on';

% actual optimization
sys_id = nlgreyest(data_merged, sys0, opt);

sys_id_vec = getpvec(sys_id);

% result visualization and saving
fit_vec = [];
[~, fit, ~] = compare(data1, sys_id);
fit_vec = [fit_vec; fit];
[~, fit, ~] = compare(data2, sys_id);
fit_vec = [fit_vec; fit];
[~, fit, ~] = compare(data3, sys_id);
fit_vec = [fit_vec; fit];
[~, fit, ~] = compare(data4, sys_id);
fit_vec = [fit_vec; fit];
disp('the worst fit to the data is: ')
disp(min(fit_vec))
disp(['M = ', num2str(sys_id_vec(1))])
disp(['m = ', num2str(sys_id_vec(2))])
disp(['b = ', num2str(sys_id_vec(3))])
disp(['c = ', num2str(sys_id_vec(4))])
param.M = sys_id_vec(1);
param.m = sys_id_vec(2);
param.b = sys_id_vec(3);
param.c = sys_id_vec(4);

compare(data1, sys_id)
