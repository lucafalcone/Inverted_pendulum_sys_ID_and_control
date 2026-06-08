u_lqi_raw = data_u_lqi(:);
u_prbs_raw = data_u_prbs(:);
y_raw = data_y_meas;

% Combine for total signal
u_total = u_lqi_raw + u_prbs_raw;

z = iddata(y_raw, u_total, 0.01);
z.InputName = {'Motor Voltage'};
z.OutputName = {'Cart Position', 'Pendulum Angle'};

z = detrend(z, 0);

% Set up sysID
sys_init = idss(A, B, C, zeros(2, 1));

% Lock rows so Matlab cant break them
sys_init.Structure.A.Free(1, :) = false;
sys_init.Structure.A.Free(3, :) = false;

sys_init.Structure.C.Free = false;
sys_init.Structure.D.Free = false;

sys_init.K = lqr(A', C', eye(4), eye(2))';
sys_init.Structure.K.Free = true;


fprintf('Running Sys ID...\n')

opt = ssestOptions('InitialState', 'zero');
sys_est = ssest(z, sys_init);
sys_est_discrete = d2d(sys_est, 0.01);

A_new = sys_est_discrete.A;
B_new = sys_est_discrete.B;

figure('Name', 'Sys ID fit verification')
compare(z, sys_est)

