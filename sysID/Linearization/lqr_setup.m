h = 0.01;
T = 60;

Q = diag([10 1 100 1]);
R = 1;
K = lqr(lin_ss, Q, R);

pause(1);
disp('done')


%% check closed loop stability
x_vec = [x v theta omega];
stable = [];
for i = 1:length(outu)
    Acl = lin_ss.A + lin_ss.B*K;
    stable = [stable; prod(eig(Acl)<0)==1]; %#ok<AGROW> 
end
plot(stable, 'x')