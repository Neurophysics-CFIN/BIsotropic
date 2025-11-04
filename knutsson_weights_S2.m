function w = knutsson_weights_S2(n, Lw, Vdeg)
% Knutsson weights restricted to even spherical-harmonic degrees.
% n : N x 3 unit vectors
% Lw: max SH degree used in the fit (only even ℓ up to Lw are included)
%
% Returns real weights w (N x 1), normalized to sum to 1.

n  = util.normalize_rows(n);
th = acos(max(-1, min(1, n(:,3))));
ph = atan2(n(:,2), n(:,1));
N  = size(n,1);

Lw = 2*floor(Lw/2); % force even
rows = {};
if nargin < 3
    k = 0:2:Lw;
    Vdeg = 1./(1+(k/2).^4);%1./(1+k.^2/36);
end
vdeg = [];
% Stack SH blocks for ℓ = 0,2,4,...,Lw   (all m each)
for l = 0:2:Lw
    Yl = util.sph_harm(l, th, ph);      % (2l+1) x N, m=-l..l (complex)
    rows{end+1} = Yl; %#ok<AGROW>
    % Per-degree weight
    v_l = Vdeg(l/2+1)/(2*l+1); %distribute weight over 2l + 1 rows
    vdeg = [vdeg; repmat(v_l, 2*l+1, 1)]; %#ok<AGROW>
end
B = vertcat(rows{:});     % M x N, M = sum_{even ℓ≤Lw}(2ℓ+1)
V = diag(vdeg);           % M x M
% Target s0: only the (ℓ=0,m=0) coefficient is nonzero.
% With Y_00 = 1/sqrt(4π), choose s0(1)=√(4π) so sum(w)≈1.
M = size(B,1);
s0 = zeros(M,1);
s0(1) = sqrt(4*pi);

% Solve (B^H V B) w = B^H V s0
A   = B' * V * B;
% Ridge for stability
lam = 1e-6 * trace(A) / max(N,1);
A = A + lam*eye(N);
rhs = B' * V * s0;

w = real(A \ rhs);

% Optional tiny negative cleanup from numerics; normalize to 1
w = max(w, 0);
s = sum(w);
if s > 0
    w = w / s;
else
    w = ones(N,1)/N;
end
end