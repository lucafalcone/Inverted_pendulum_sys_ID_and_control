function predmod = predmodgen(LTI, dim)
nx = dim.nx;
nu = dim.nu;
N  = dim.N;

%% Prediction matrix from initial state
% T = [A; A^2; ...; A^N]
T = zeros(N*nx, nx);
for k = 1:N
    T((k-1)*nx+1 : k*nx, :) = LTI.A^k;
end

%% Prediction matrix from input sequence
% S is lower-triangular block Toeplitz:
%   S(i,j) = A^(i-j) * B  for i >= j,  0 otherwise
S = zeros(N*nx, N*nu);
for i = 1:N
    for j = 1:i
        S((i-1)*nx+1 : i*nx, (j-1)*nu+1 : j*nu) = LTI.A^(i-j) * LTI.B;
    end
end

predmod.T = T;
predmod.S = S;

end
