proj_root = fileparts(fileparts(mfilename('fullpath')));
load(fullfile(proj_root, 'LQI', 'results', 'lqr_kalman.mat'));

T = 30;
N_samples = T / h;
t_prbs = (0:N_samples-1)' * h;

% Generate the PRBS signal
prbs_data = idinput(N_samples, 'PRBS', [0 0.2], [-1 1]);

prbs_ts = timeseries(prbs_data, t_prbs);

%%
mname = 'inverted_pendulum_LQR';
play_run(mname)