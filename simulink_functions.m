%% PREPARE WORKSPACE
clear
clc

exptype = 3;
h = 0.01;
T = 120;

t = (0:h:T)';

switch exptype
    case 1
        inputtype = '05_1_Doub';
        A = 0.5;
        T_doublet = 1;
        u = @(t) A .* (t >= 0 & t < T_doublet) + ...
            -A .* (t >= T_doublet & t < 2*T_doublet);
        u = timeseries(u(t),t);
    case 2
        inputtype = '1_05_Doub';
        A = 1;
        T_doublet = 0.5;
        u = @(t) A .* (t >= 0 & t < T_doublet) + ...
            -A .* (t >= T_doublet & t < 2*T_doublet);
        u = timeseries(u(t),t);
    case 3
        inputtype = '100_005_Doub';
        A = 100;
        T_doublet = 0.05;
        u = @(t) A .* (t >= 0 & t < T_doublet) + ...
            -A .* (t >= T_doublet & t < 2*T_doublet);
        u = timeseries(u(t),t);
    case 4
        inputtype = 'step';
        A = -0.5;
        u = @(t) A .* (t >= 1 & t < 10);
        u = timeseries(u(t),t);
end



disp('...')
pause(1)
disp('done loading')

%% SAVE OUTPUT DATA 
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
filename = sprintf('%s_%s.mat', inputtype, timestamp);
save(fullfile(pwd, 'real_data', filename), 'outu', 'x', 'theta')

disp('saved data')