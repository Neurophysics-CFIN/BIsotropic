function val = electroSO3_cost_fct(x, N,addpen)
%   ELECTROSO3_COST_FCT  Electrostatic-style repulsion energy for SO(3) samples (quaternions).
%
%   val = ELECTROSO3_COST_FCT(x, N)
%   val = ELECTROSO3_COST_FCT(x, N, addpen)
%
%   Computes a repulsion (Coulomb-like) cost for a set of N rotations
%   represented as unit quaternions (on the 3-sphere, S^3). The cost is based on pairwise SO(3)
%   geodesic distances derived from quaternion inner products with antipodal
%   identification (q ≡ -q). Optionally, additional penalties can be added
%   to discourage close pairs and encourage distance uniformity.
%
%   INPUTS
%     x      : (4N-by-1) vector containing N quaternions stacked as rows
%              [t x y z]. Reshaped to Q = reshape(x,[N,4]).
%     N      : Number of samples/quaternions.
%     addpen : (optional) Logical flag enabling additional penalty terms.
%              If omitted or false, only the base mean inverse-distance
%              energy is used.
%
%   OUTPUT
%     val : Scalar cost value.
%
%   DETAILS
%     - Q is row-normalized to unit quaternions and canonicalized so that
%       q and -q map to a unique representative (SO(3) folding).
%     - Pairwise similarities are computed as S = |Q*Q'|, clamped to [0,1].
%     - Pairwise normalized geodesic distances on SO(3) are then
%           Phi = (2/pi) * acos(S),
%       which lies in [0,1].
%     - Self-distances are ignored and only unordered pairs (i<j) are used.
%     - Base cost (default):
%           val = mean( 1 ./ (Phi + eps) ).
%
%   OPTIONAL PENALTIES (when addpen is true)
%     - Close-pair sharpening: uses exponent p=0.5 in
%           mean( 1 ./ (Phi + eps).^p )
%       (with eps = 1e-12 for numerical safety).
%     - Uniformity penalty: adds lambdaVar * var(Phi) over off-diagonal
%       distances, with lambdaVar = 100.
%
%   NOTES
%     - The distance normalization (2/pi) maps SO(3) rotation angle distance
%       (in [0,pi]) to [0,1].
%     - The function assumes util.normalize_rows and util.canon_antipodes
%       are available on the MATLAB path.
%
%   EXAMPLE
%     N = 300;
%     x = randn(4*N,1);
%     val0 = electroSO3_cost_fct(x, N);
%     val1 = electroSO3_cost_fct(x, N, true);


Q = reshape(x, [N,4]);

Q = util.normalize_rows(Q);
Q = util.canon_antipodes(Q);

S = abs(Q*Q.');                         % |dot| in [0,1] (up to fp error)
S = min(1, max(0, S));                  % clamp numerically
Phi = (2/pi) * acos(S);                 % normalized geodesic distances
n = size(Q,1);
Phi(1:n+1:end) = NaN;                   % ignore self-distances
w = Phi(triu(true(n),1));               % upper triangle


% additional penalizing terms
% Hyperparameters
if nargin == 4 && addpen
    p           = .5;         % sharper close-pair penalty 2
    lambdaVar   = 100;       % distance variance weight .2

    % Uniformity penalty (variance of distances, off-diagonal)
    varD = var(w,1);

    val = mean(1 ./ (w(~isnan(w)) + 1e-12).^p);
    val = val + lambdaVar*varD;

else
    val = mean(1 ./ (w(~isnan(w)) + 1e-12));
end

end

