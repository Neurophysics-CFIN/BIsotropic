function [t,x,y,z] = hopf_kronecker(N)
%HOPF_KRONECKER  Hopf–Kronecker lattice on S^3 for deterministic SO(3) sampling.
%
%   T = HOPF_KRONECKER(N)
%   [t,x,y,z] = HOPF_KRONECKER(N)
%
%   Generates N deterministic unit quaternions on S^3 using a Hopf–Kronecker
%   (golden-ratio / irrational rotation) lattice, then folds antipodes
%   (q ≡ -q) to obtain a canonical representative suitable for sampling SO(3).
%
%   INPUT
%     N : Positive integer number of samples.
%
%   OUTPUT
%     T : (N-by-4) array of unit quaternions (versors), one per row, ordered as
%         [t x y z]. Returned when called with one output.
%
%     t,x,y,z : Each is (N-by-1). The quaternion components returned as
%               separate vectors when called with four outputs. 
% 
%               t scalar part
%                  
%
%   DETAILS
%     - Uses two incommensurate angles based on a = (sqrt(5)-1)/2 and
%       b = sqrt(2)-1 to avoid periodicity and produce well-distributed
%       phases on the Hopf torus fibers.
%     - Uses u in (-1,1) and psi = 0.5*acos(u) to spread samples in the
%       Hopf coordinate "latitude" direction.
%     - Applies a lexicographic sign rule to choose a canonical element
%       from each antipodal pair {q, -q}, matching the SO(3) identification.
%     - Normalizes as a safety step to ensure unit quaternions.
%
%   EXAMPLE
%     % Get quaternions as an N-by-4 matrix
%     Q = hopf_kronecker(1000);
%
%     % Or get components separately
%     [t,x,y,z] = hopf_kronecker(1000);

i = (0:N-1)'; a = (sqrt(5)-1)/2; b = sqrt(2)-1;
phi1 = 2*pi*mod(a*i,1);  phi2 = 2*pi*mod(b*i,1);
u = -1 + 2*(i+0.5)/N;    psi = 0.5*acos(u);
C = cos(psi); S = sin(psi);
Q = [C.*cos(phi1), C.*sin(phi1), S.*cos(phi2), S.*sin(phi2)]; % quats

% Canonicalize antipodes for SO(3)
neg = Q(:,1)<0 | (Q(:,1)==0 & (Q(:,2)<0 | (Q(:,2)==0 & (Q(:,3)<0 | (Q(:,3)==0 & Q(:,4)<0)))));
Q(neg,:) = -Q(neg,:);
Q = Q ./ vecnorm(Q,2,2); % safety normalize
if nargout ==1
    t = Q;
else
    t = Q(:,1); x = Q(:,2); y = Q(:,3); z = Q(:,4);
end
end