function [A, B, C, sv] = okid_era(u, y, nx, p)
% OKID_ERA  Black-box subspace identification via OKID + ERA.
%   Identifies a discrete-time LTI model from closed-loop I/O data.
%   Assumes zero direct feedthrough (D = 0), which holds for the inverted
%   pendulum (motor command does not instantaneously appear in position).
%
%   Inputs
%     u   : N×m  input data  (e.g. motor command 'outu')
%     y   : N×l  output data (e.g. [x, theta])
%     nx  : desired model order
%     p   : number of observer Markov parameters (default: min(floor(N/4),60))
%
%   Outputs
%     A, B, C  : discrete-time state-space matrices (at the sample rate of u/y)
%     sv       : singular values from block-Hankel SVD (plot to verify order)
%
%   Algorithm
%     OKID  – fits an ARX-like observer model to the closed-loop data.
%             Including past outputs in the regressor makes the regression
%             stable even when the open-loop plant is unstable.
%     ERA   – builds a block-Hankel matrix from the recovered system impulse
%             response, SVD-truncates to nx states, and extracts A, B, C.

[N, m] = size(u);
l = size(y, 2);
if nargin < 4
    p = min(floor(N / 5), 100);
end

% ── OKID regression ───────────────────────────────────────────────────────
% Observer model (D = 0):
%   y(k) = Σ_{i=1}^{p} M_i*u(k-i) + Σ_{i=1}^{p} L_i*y(k-i)
%
% Regressor at step k (1-indexed, k > p):
%   v(k) = [u(k-1); y(k-1); u(k-2); y(k-2); ...; u(k-p); y(k-p)]
%   length = p*(m+l)

nreg = p * (m + l);
Nv   = N - p;

V     = zeros(nreg, Nv);
Y_mat = zeros(l, Nv);

for k = 1:Nv
    kk = k + p;
    v  = zeros(nreg, 1);
    for i = 1:p
        v((i-1)*(m+l)+1   : (i-1)*(m+l)+m)   = u(kk-i, :)';
        v((i-1)*(m+l)+m+1 : i*(m+l))          = y(kk-i, :)';
    end
    V(:, k)     = v;
    Y_mat(:, k) = y(kk, :)';
end

% Least-squares solve with Tikhonov regularisation
lambda = 1e-6 * trace(V * V') / nreg;
Theta  = Y_mat * V' / (V * V' + lambda * eye(nreg));   % l × nreg

% ── Recover system Markov parameters ─────────────────────────────────────
% Recursion (D = 0, so h_0 = 0):
%   h_k = M_k + Σ_{j=1}^{k-1} L_j * h_{k-j}
%
% M_k and L_k are the input and output blocks of Theta at lag k.

h = zeros(l, m, p);   % h(:,:,k) = system Markov parameter at lag k

for k = 1:p
    off     = (k-1)*(m+l);
    M_k     = Theta(:, off+1   : off+m);
    h(:,:,k) = M_k;
    for j = 1:k-1
        off_j    = (j-1)*(m+l);
        Lj       = Theta(:, off_j+m+1 : off_j+m+l);
        h(:,:,k) = h(:,:,k) + Lj * h(:,:,k-j);
    end
end

% ── ERA: block-Hankel SVD → state-space ───────────────────────────────────
nr = floor(p * 0.60);   % block rows    (nr*l ≥ nx required)
nc = p - nr - 1;        % block columns (nc*m ≥ nx required, H1 uses lag nr+nc+1 ≤ p)

if nr * l < nx || nc * m < nx
    error('p=%d too small for nx=%d. Increase p.', p, nx);
end

H0 = zeros(nr*l, nc*m);
H1 = zeros(nr*l, nc*m);
for r = 1:nr
    for c = 1:nc
        H0((r-1)*l+1:r*l, (c-1)*m+1:c*m) = h(:,:,r+c-1);
        H1((r-1)*l+1:r*l, (c-1)*m+1:c*m) = h(:,:,r+c);
    end
end

[U_s, S_s, V_s] = svd(H0);
sv = diag(S_s);

if nx > length(sv)
    error('nx=%d exceeds available singular values (%d). Reduce nx or increase p.', ...
          nx, length(sv));
end

sq  = sqrt(sv(1:nx));
isq = 1 ./ sq;
U1  = U_s(:, 1:nx);
V1  = V_s(:, 1:nx);

O = U1 * diag(sq);      % observability matrix:   nr*l × nx
R = diag(sq) * V1';     % controllability matrix: nx × nc*m

C = O(1:l, :);          % l  × nx
B = R(:, 1:m);          % nx × m
A = diag(isq) * (U1' * H1 * V1) * diag(isq);   % nx × nx
end
