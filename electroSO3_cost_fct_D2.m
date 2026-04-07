function val = electroSO3_cost_fct_D2(x, N)
%   ELECTROSO3_COST_FCT_D2  D2-symmetric variant of ELECTROSO3_COST_FCT via right Q8 quotient.
%
%   val = ELECTROSO3_COST_FCT_D2(x, N)
%
%   Same base electrostatic SO(3) repulsion cost as ELECTROSO3_COST_FCT
%   (mean inverse normalized geodesic distance from quaternion dot products),
%   but modified to enforce right dihedral (D2) symmetry by quotienting the
%   sample set under the right action of the quaternion group Q8.
%
%   KEY DIFFERENCE FROM ELECTROSO3_COST_FCT
%     - Instead of using S = |Q*Q'| directly, this routine accounts for the
%       symmetry by right-multiplying the second factor by each K in Q8 and
%       taking the elementwise maximum similarity:
%           S(j,k) = max_{K in Q8} | < q_j , q_k ⊗ K > |.
%       This effectively identifies orientations that differ by a right D2
%       symmetry operation (covered by Q8 in the quaternion lift).
%
%   Everything else (normalization, antipodal canonicalization, conversion
%   to normalized geodesic distance Phi = (2/pi)*acos(S), dropping self-terms,
%   and the mean of 1/(Phi+eps) over unordered pairs) matches
%   ELECTROSO3_COST_FCT.
%
%   SEE ALSO
%     ELECTROSO3_COST_FCT

Q = reshape(x, [N,4]);
%test
% Q = [1 0 0 0;Q];
Q = util.normalize_rows(Q);
Q = util.canon_antipodes(Q);

% --- Quotient by right action of Q8 (covers D2 symmetry in SO(3)) ---
% Q8 elements in convention [t x y z]  (t = scalar part)
K = [ 1  0  0  0;   %  +1
     -1  0  0  0;   %  -1  (redundant with abs(.), but cheap)
      0  1  0  0;   %  +i
      0 -1  0  0;   %  -i
      0  0  1  0;   %  +j
      0  0 -1  0;   %  -j
      0  0  0  1;   %  +k
      0  0  0 -1];  %  -k

n = size(Q,1);
S = zeros(n,n);

for kk = 1:size(K,1)
    Qk = qmul_right(Q, K(kk,:));   % Qk(j,:) = Q(j,:) ⊗ K(kk,:)
    S  = max(S, abs(Q * Qk.'));    % elementwise max over kk
end
S = min(1, max(0, S));             % clamp numerically

Phi = (2/pi) * acos(S);                 % normalized geodesic distances
n = size(Q,1);
Phi(1:n+1:end) = NaN;                   % ignore self-distances
w = Phi(triu(true(n),1));               % upper triangle
val = mean(1 ./ (w(~isnan(w)) + 1e-12));
end

function Qout = qmul_right(Q, r)
%Qout = Q ⊗ r, with quaternions in [t x y z]
% Q: n×4, r: 1×4
t = Q(:,1); x = Q(:,2); y = Q(:,3); z = Q(:,4);
rt = r(1);  rx = r(2);  ry = r(3);  rz = r(4);

Qout = zeros(size(Q),'like',Q);
Qout(:,1) = t*rt - x*rx - y*ry - z*rz;
Qout(:,2) = t*rx + x*rt + y*rz - z*ry;
Qout(:,3) = t*ry - x*rz + y*rt + z*rx;
Qout(:,4) = t*rz + x*ry - y*rx + z*rt;
end
