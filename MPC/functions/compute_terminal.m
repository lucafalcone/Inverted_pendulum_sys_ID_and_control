function terminal = compute_terminal(LTI, param)
%% Solve DARE for terminal cost P and LQR gain K
[P, ~, Kd] = dare(LTI.A, LTI.B, param.Q, param.R);
K = -Kd;  % dare returns K such that u = -Kd*x; we use u = K*x with K = -Kd

terminal.P = P;
terminal.K = K;

%% Verify closed-loop stability
A_cl = LTI.A + LTI.B * K;
eig_cl = eig(A_cl);
fprintf('Terminal controller eigenvalues:\n');
for i = 1:length(eig_cl)
    fprintf('  lambda_%d = %.4f + %.4fi  (|lambda| = %.4f)\n', ...
        i, real(eig_cl(i)), imag(eig_cl(i)), abs(eig_cl(i)));
end
if all(abs(eig_cl) < 1)
    fprintf('  => Closed-loop system is asymptotically stable.\n\n');
else
    warning('Closed-loop system is NOT stable!');
end

%% Compute terminal invariant set X_f
% The terminal set is the maximal positively invariant set for the
% closed-loop system x(k+1) = A_cl * x(k) subject to:
%   (a) Input constraints:  |K*x| <= u_max  (component-wise)
%   (b) State constraints:  approach cone constraints
%
% We compute X_f iteratively using the algorithm from Rawlings & Mayne.

% Constraint matrices for input: |K*x| <= u_max
% F_u = [K; -K];
% g_u = limits.u_max * ones(2*limits.nu, 1);
% 
% % Constraint matrices for approach cone
% tg = tan(param.gamma_los);
% F_x = [ tg,  1, 0, 0;
%         tg, -1, 0, 0];
% g_x = [0; 0];
% 
% % Combined constraint: F * x <= g
% F = [F_u; F_x];
% g = [g_u; g_x];
% 
% Iterative computation of maximal positively invariant set
% Iterative computation of the maximal positively invariant set.
% At each iteration, we propagate constraints through the closed-loop
% % dynamics and check convergence by comparing the constraint set norm.
% F_set = F;
% g_set = g;
% 
% max_iter = 100;
% converged = false;
% 
% for iter = 1:max_iter
%     % Add constraints for next step: F * A_cl^iter * x <= g
%     F_new = F * A_cl^iter;
% 
%     % Check convergence: if new constraints are negligible (rows decay
%     % because A_cl is stable), we have converged.
%     if norm(F_new, 'fro') < 1e-10
%         converged = true;
%         fprintf('Terminal set converged after %d iterations (row decay).\n', iter);
%         break;
%     end
% 
%     F_set = [F_set; F_new];
%     g_set = [g_set; g];
% end
% 
% if ~converged
%     fprintf('Terminal set computation stopped after %d iterations.\n', max_iter);
% end
% 
% fprintf('Terminal set has %d inequality constraints.\n', size(F_set, 1));
% 
% terminal.Xf_A = F_set;
% terminal.Xf_b = g_set;

end
