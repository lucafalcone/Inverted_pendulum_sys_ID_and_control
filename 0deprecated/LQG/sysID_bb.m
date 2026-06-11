%% Black-box closed-loop system identification from LQG chirp experiment
% Fits a 4th-order discrete state-space model using n4sid + ssest on the
% data collected by main_LQG.m (LQG stabilises upright, chirp perturbs it).
%
% Plant input:  F_total = F_cmd + delta_F  [N]  (force after dead-zone block)
% Plant output: [x, theta]  [m, rad]
%
% Run after main_LQG.m has completed a hardware experiment.

clear; clc;
h  = 0.01;  % sample time [s]
nx = 4;     % model order: [x, v, theta, omega]

proj_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(proj_root, 'global'));

results_dir = fullfile(proj_root, 'LQG', 'results');
timestamp   = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

%% 1) Load latest hardware run
run_dirs = dir(fullfile(results_dir, '2*'));
if isempty(run_dirs)
    error('No results in LQG/results. Run main_LQG.m first.');
end
[~, idx] = max([run_dirs.datenum]);
run_dir  = fullfile(results_dir, run_dirs(idx).name);
files    = dir(fullfile(run_dir, 'run_hw_*.mat'));
if isempty(files)
    error('No run_hw_*.mat found in %s', run_dir);
end
S = load(fullfile(run_dir, files(1).name));
fprintf('Loaded: %s\n', fullfile(run_dir, files(1).name));

%% 2) Extract and align signals
% Drop the first second: controller hasn't fully settled and chirp starts
% at near-zero amplitude, contributing little information.
skip    = round(1.0 / h);
F_total = S.F_total(skip+1:end);
x       = S.x(skip+1:end);
theta   = S.theta(skip+1:end);

% Trim to equal length (ToWorkspace blocks can differ by 1 sample)
N = min([length(F_total), length(x), length(theta)]);
F_total = F_total(1:N);
x       = x(1:N);
theta   = theta(1:N);

fprintf('Using %d samples (%.1f s) after skip\n', N, N*h);

%% 3) Split: 70 % estimation / 30 % validation
N_est = floor(0.70 * N);
N_val = N - N_est;

make_dat = @(u, y) iddata(y, u, h, ...
    'OutputName', {'x'; 'theta'}, 'InputName', {'F'}, ...
    'OutputUnit', {'m'; 'rad'},   'InputUnit',  {'N'});

dat_est = detrend(make_dat(F_total(1:N_est),           [x(1:N_est),     theta(1:N_est)]),     0);
dat_val = detrend(make_dat(F_total(N_est+1:end), [x(N_est+1:end), theta(N_est+1:end)]), 0);

%% 4) n4sid — subspace initialisation
fprintf('\n--- n4sid (order %d) ---\n', nx);
opt_n4 = n4sidOptions('Focus', 'prediction', 'EnforceStability', false);
sys_n4 = n4sid(dat_est, nx, opt_n4);

fprintf('n4sid poles (discrete):\n');
disp(eig(sys_n4.A));

%% 5) ssest — prediction-error refinement starting from n4sid
fprintf('--- ssest refinement ---\n');
opt_ss = ssestOptions('InitialState', 'estimate', 'EnforceStability', false);
sys_bb = ssest(dat_est, sys_n4, opt_ss);

fprintf('ssest poles (discrete):\n');
p_bb = eig(sys_bb.A);
disp(p_bb);
if any(abs(p_bb) > 1)
    fprintf('  -> Unstable pole confirmed (expected)\n');
end

%% 6) Validation fit on held-out data
[~, fit] = compare(dat_val, sys_bb);
fprintf('\nValidation fit (30%% held-out):\n');
fprintf('  x     : %.1f%%\n', fit(1));
fprintf('  theta : %.1f%%\n', fit(2));

%% 7) Compare against grey-box linearisation
p_gb = load(fullfile(proj_root, 'global', 'param_64_175.mat')).param;
p_gb.g = 9.81;
f_eom  = @(z, F) eom_force(z, F, p_gb.M, p_gb.m, p_gb.b, p_gb.c, p_gb.l);
z0 = zeros(4,1); eps_j = 1e-6;
A_c = zeros(4,4);
for i = 1:4
    dz = zeros(4,1); dz(i) = eps_j;
    A_c(:,i) = (f_eom(z0+dz,0) - f_eom(z0-dz,0)) / (2*eps_j);
end
B_c     = (f_eom(z0, eps_j) - f_eom(z0, -eps_j)) / (2*eps_j);
sys_phys = c2d(ss(A_c, B_c, [1 0 0 0; 0 0 1 0], zeros(2,1)), h);

fprintf('\nGrey-box poles (discrete):\n');
disp(eig(sys_phys.A));

%% 8) Plots
% Validation: compare predicted vs measured on held-out data
fig_val = figure('Name', 'Validation: black-box vs measured');
compare(dat_val, sys_bb);
title('Black-box model vs measured — validation set (30%)');
saveas(fig_val, fullfile(run_dir, ['bb_validation_' timestamp '.png']));

% Bode comparison
w = logspace(-1, 2, 300);
fig_bode = figure('Name', 'Bode: black-box vs grey-box');
bodeplot(sys_bb, sys_phys, w);
legend('Black-box (n4sid+ssest)', 'Grey-box linearisation', 'Location', 'best');
title('Bode: identified vs physics model');
grid on;
saveas(fig_bode, fullfile(run_dir, ['bode_comparison_' timestamp '.png']));

% Pole comparison
fig_pz = figure('Name', 'Poles: black-box vs grey-box');
pzplot(sys_bb, sys_phys);
legend('Black-box', 'Grey-box');
title('Poles and zeros — identified vs physics model');
grid on;
saveas(fig_pz, fullfile(run_dir, ['pzmap_' timestamp '.png']));

%% 9) Save
out.h        = h;
out.sys_bb   = sys_bb;
out.A        = sys_bb.A;  out.B = sys_bb.B;  out.C = sys_bb.C;
out.sys_phys = sys_phys;
out.fit      = fit;
out.dat_est  = dat_est;
out.dat_val  = dat_val;

save(fullfile(run_dir,    ['bb_model_' timestamp '.mat']), '-struct', 'out');
save(fullfile(results_dir, 'bb_model_latest.mat'),         '-struct', 'out');
fprintf('\nSaved to LQG/results/bb_model_latest.mat\n');
