function val = GFOSO3_cost_fct(x, N, Lmax, bandWeight)
%   GFOSO3_COST_FCT  Band-weighted SO(3) sampling energy using quaternion dot products.
%
%   val = GFOSO3_COST_FCT(x, N, Lmax, bandWeight)
%
%   Computes a scalar objective (energy/cost) for an SO(3) point set encoded
%   as quaternions. The input vector x is reshaped into N quaternions, each
%   row is normalized to unit length, and a band-limited pairwise energy is
%   accumulated using Chebyshev polynomials of the second kind evaluated on
%   the pairwise quaternion inner-product matrix, corresponding to SO(3) 
%   group characters.
%
%   INPUTS
%     x          : (4N-by-1) vector representing N quaternions stacked in
%                  row-major order. Reshaped as Q = reshape(x,[N,4]), where
%                  each row is [t x y z], t scalar part.
%     N          : Number of quaternions/samples.
%     Lmax       : Maximum bandlimit. Bands ℓ = 2...Lmax are accumulated.
%     bandWeight : Band weights (V_{ll}^2). The code uses bandWeight(ell+1) for
%                  ell = 2..Lmax, so bandWeight must have length >= Lmax+1.
%                  (Entries for ℓ=0 and ℓ=1 may exist and are intentionally
%                  unused by this routine.)
%
%   OUTPUT
%     val : Scalar cost value.
%
%   DETAILS
%     - Normalization: Q is row-normalized to enforce unit quaternions.
%     - Pairwise matrix: S = Q*Q' is clamped to (-1,1) for numerical
%       stability and its diagonal is set to 1.
%     - Recurrence: Chebyshev polynomials of the second kind U_k are
%       advanced entrywise via
%           U_{k+1}(S) = 2*S .* U_k(S) - U_{k-1}(S),
%       with U_0 = 1 and U_1 = 2*S.
%     - Accumulation: Only even k are used, with ell = k/2. For each ell,
%       the routine sums all off-diagonal entries of U_{2*ell}(S) and adds
%       bandWeight(ell+1) times this sum to val. Diagonal (self-pair) terms
%       are intentionally excluded as they do depend on Q.
%
%   NOTES
%     - Antipodal folding (q ≡ -q) is intentionally not applied here, but
%     should be applied in calling routine, e.g. isotropic_SO3.m
%       (canonicalization is shown but commented out in the code).
%     - The current implementation uses the signed dot products S (not |S|), 
%       as U_{2k} is even.
%
%   EXAMPLE
%     N = 500; Lmax = 12;
%     bandWeight = util.sobolev(8,7,Lmax);
%     x = randn(4*N,1);
%     val = GFOSO3_cost_fct(x, N, Lmax, bandWeight);

% x: 4N×1 vector → N×4 quaternions (rows). We normalize rows and canonicalize.

Q = reshape(x, [N,4]);
Q = util.normalize_rows(Q);              % enforce unit quaternions
% Q = util.canon_antipodes(Q);            % q ≡ -q
nR = size(Q,1);

% Pairwise |dot| using sqrt(s^2) (no abs; stable)
S = Q*Q.';                 % nR×nR
% S = sqrt(S.^2);            % entrywise |.| in [0,1]

S = max(min(S, 1-1e-12), -1+1e-12);
S(1:nR+1:end) = 1;
U_prev = ones(nR);     % U0
U_curr = 2*S;          % U1

% advance to U2 (k=2), no contribution (ℓ=1)
U_next = 2*S .* U_curr - U_prev; 
U_prev = U_curr; U_curr = U_next;

val = 0;
for k = 3:2*Lmax
    U_next = 2*S .* U_curr - U_prev;
    if mod(k,2)==0
        ell = k/2;                    % ell = 2..Lmax
        sum_all = sum(U_next(:)) - sum(diag(U_next));   % subtracting constant diagonal
        val = val + bandWeight(ell+1) * sum_all;
    end
    U_prev = U_curr; U_curr = U_next;
end

end