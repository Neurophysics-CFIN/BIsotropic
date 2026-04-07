function Q = canon_antipodes(Q)
%CANON_ANTIPODES  Choose a canonical representative for antipodal points on
%2- or 3-sphere.
%
%   Q = CANON_ANTIPODES(Q)
%
%   Selects a unique, canonical representative from each antipodal pair
%   q ≡ -q, as commonly required when working with SO(3) representations.
%   The choice is made using a lexicographic tie-break rule so that results
%   are deterministic and consistent.
%
%   INPUT
%     Q : An N-by-4 or N-by-3 array.
%         - If N-by-4, rows are unit quaternions [t x y z] (SO(3)).
%         - If N-by-3, rows are 3-vectors with antipodal symmetry (S^2).
%
%   OUTPUT
%     Q : The modified array where each row has been flipped, if necessary,
%         so that it is the canonical representative of {q, -q}.
%
%   NOTES
%     - For quaternions, the first nonzero component is enforced to be
%       positive (lexicographic ordering).
%     - This is useful for SO(3) sampling, clustering, and comparisons,
%       where q and -q represent the same rotation.
%
%   EXAMPLE
%     % Canonicalize a set of random quaternions
%     Q = randn(10,4);
%     Q = Q ./ vecnorm(Q,2,2);
%     Qc = canon_antipodes(Q);

% Pick canonical representative of q ≡ -q (lexicographic tie-break)
if size(Q,2) == 4
    neg = Q(:,1) < 0 | (Q(:,1)==0 & (Q(:,2) < 0 | (Q(:,2)==0 & (Q(:,3) < 0 | (Q(:,3)==0 & Q(:,4) < 0)))));
    Q(neg,:) = -Q(neg,:);
elseif size(Q,2) ==3
    neg = Q(:,3) < 0 | (Q(:,3)==0 & (Q(:,2) < 0 | (Q(:,2)==0 & (Q(:,1) < 0 ))));
    Q(neg,:) = -Q(neg,:);
else
    error('Input matrix Q must have 3 or 4 columns.');
end
