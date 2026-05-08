clear
clc

% creating id data object
load real_data\fixedCart_20260506_155243.mat
data1 = iddata(theta(1:500), [], (t(2) - t(1)));
data2 = iddata(theta(1000:1500), [], (t(2) - t(1)));
load real_data\fixedCart_20260506_175042.mat
data3 = iddata(theta(1:500), [], (t(2) - t(1)));
data4 = iddata(theta(1000:1500), [], (t(2) - t(1)));
data_merged = merge(data1, data2, data3, data4);


% creating the object for parameter estimation
order = [1 0 2];
params0 = [0.05, 0.3];     % b/m and l
x0 = [theta(1) 0]';

sys0 = idnlgrey('fixed_cart_EOM', order, params0, x0, 0);

sys0.Parameters(1).Minimum = 0;
sys0.Parameters(2).Minimum = 0;
sys0.InitialStates(1).Fixed = false;
sys0.InitialStates(2).Fixed = false;

opt = nlgreyestOptions;
opt.SearchMethod = 'auto';
opt.SearchOptions.MaxIterations = 30;
opt.Display = 'on';

% actual optimization
sys_id = nlgreyest(data_merged, sys0, opt);

sys_id_vec = getpvec(sys_id);

% result visualization
fit_vec = [];
[~, fit, ~] = compare(data1, sys_id);
fit_vec = [fit_vec; fit];
[~, fit, ~] = compare(data1, sys_id);
fit_vec = [fit_vec; fit];
[~, fit, ~] = compare(data1, sys_id);
fit_vec = [fit_vec; fit];
[~, fit, ~] = compare(data1, sys_id);
fit_vec = [fit_vec; fit];
disp('the worst fit to the data is: ')
disp(min(fit_vec))
disp('values for b/m and for l are:')
disp(sys_id_vec)

compare(data1, sys_id)
