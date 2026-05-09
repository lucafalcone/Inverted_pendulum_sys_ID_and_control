%% PREPARE WORKSPACE
clear
clc

exptype = 11;
h = 0.01;
T = 20;

t = (0:h:T)';

switch exptype
    case 0
        inputtype = 'step';
        A = 0;
        duration = 2;
        u = @(t) A .* (t >= 0 & t < duration);
    case 1
        inputtype = '05_1_Doub';
        A = 0.5;
        T_doublet = 1;
        u = @(t) A .* (t >= 0 & t < T_doublet) + ...
            -A .* (t >= T_doublet & t < 2*T_doublet);
    case 2
        inputtype = '1_05_Doub';
        A = 1;
        T_doublet = 0.5;
        u = @(t) A .* (t >= 0 & t < T_doublet) + ...
            -A .* (t >= T_doublet & t < 2*T_doublet);
    case 3
        inputtype = '100_005_Doub';
        A = 100;
        T_doublet = 0.05;
        u = @(t) A .* (t >= 0 & t < T_doublet) + ...
            -A .* (t >= T_doublet & t < 2*T_doublet);
    case 4
        inputtype = 'step';
        A = -0.5;
        u = @(t) A .* (t >= 1 & t < 10);
    case 5
        inputtype = 'fixedCart';
        u = @(t) t .* 0;
    case 6
        inputtype = '1sin1';
        u = @(t) 1*sin(1*pi*t);
    case 7
        inputtype = '1sin2';
        u = @(t) 1*sin(2*pi*t);
    case 8
        inputtype = '1sin4';
        u = @(t) 1*sin(4*pi*t);
    case 9
        inputtype = 'multisin';
        u = @(t) 0.5*sin(0.7*pi*t) + 0.5*sin(1.25*pi*t) + 0.5*sin(3.24*pi*t);
    case 10 
        inputtype = 'chirp1';
        u = @(t) chirp(t,0.5,20,4,'linear');    
    case 11
        inputtype = 'chirp2';
        u = @(t) chirp(t,0.25,20,3,'linear');
end
u = timeseries(u(t),t);


disp('...')
pause(1)
disp('done loading')

%% SAVE OUTPUT DATA 
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
filename = sprintf('%s_%s.mat', inputtype, timestamp);
x = x - x(1);
save(fullfile(pwd, 'real_data', filename), 'outu', 'x', 'theta', 't')

disp('saved data')