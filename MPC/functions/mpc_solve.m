function [u_opt, fval, exitflag] = mpc_solve(x0, H, h, constraints, dim)
%% Set up QP (augmented with slack s for soft position constraint)
% z = [U (N*nu x 1); s (scalar)], penalty on s is in H(end,end)
f      = [h * x0; 0];   % zero linear cost on slack
b_rhs  = constraints.b_ineq + constraints.bx0 * x0;

%% Solve with quadprog
options = optimoptions('quadprog', ...
    'Display',   'none', ...
    'Algorithm', 'interior-point-convex');

[z_opt, fval, exitflag] = quadprog(H, f, constraints.A_ineq, b_rhs, ...
    [], [], [], [], [], options);

if exitflag ~= 1
    warning('quadprog exitflag = %d at state [%.2f, %.2f, %.2f, %.2f]', ...
        exitflag, x0(1), x0(2), x0(3), x0(4));
    if isempty(z_opt)
        z_opt = zeros(dim.N * dim.nu + 1, 1);
    end
end

u_opt = z_opt(1);

end
