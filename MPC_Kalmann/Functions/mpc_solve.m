function [u_opt, fval, exitflag] = mpc_solve(x0, H, h, constraints, dim)
%% Set up QP
f      = h * x0;
b_rhs  = constraints.b_ineq + constraints.bx0 * x0;

%% Solve with quadprog
options = optimoptions('quadprog', ...
    'Display',   'none', ...
    'Algorithm', 'interior-point-convex');

[u_opt, fval, exitflag] = quadprog(H, f, constraints.A_ineq, b_rhs, ...
    [], [], [], [], [], options);

if exitflag ~= 1
    warning('quadprog exitflag = %d at state [%.2f, %.2f, %.2f, %.2f]', ...
        exitflag, x0(1), x0(2), x0(3), x0(4));
    if isempty(u_opt)
        u_opt = zeros(dim.N * dim.nu, 1);
    end
end

end
