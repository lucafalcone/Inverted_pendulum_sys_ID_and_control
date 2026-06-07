function [H, h, const] = costgen(predmod, weight, dim)
nx = dim.nx;
N  = dim.N;

%% Build block-diagonal weight matrix
% Qbar = blkdiag(Q, Q, ..., Q, P)  size (N*nx) x (N*nx)
Qbar = blkdiag(kron(eye(N-1), weight.Q), weight.P);
Rbar = kron(eye(N), weight.R);

%% Quadratic cost matrices
% J = 0.5 * U' * H * U + (h * x0)' * U + x0' * const * x0
H = predmod.S' * Qbar * predmod.S + Rbar;
H = (H + H') / 2;  % Ensure symmetry

h = predmod.S' * Qbar * predmod.T;

const = predmod.T' * Qbar * predmod.T;

end
