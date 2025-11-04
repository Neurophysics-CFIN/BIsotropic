function [t,x,y,z] = hopf_kronecker(N)
%% ============================================================
%  Hopf–Kronecker S^3 lattice → SO(3) via antipodal folding
% =============================================================
    i = (0:N-1)'; a = (sqrt(5)-1)/2; b = sqrt(2)-1;
    phi1 = 2*pi*mod(a*i,1);  phi2 = 2*pi*mod(b*i,1);
    u = -1 + 2*(i+0.5)/N;    psi = 0.5*acos(u);
    C = cos(psi); S = sin(psi);
    Q = [C.*cos(phi1), C.*sin(phi1), S.*cos(phi2), S.*sin(phi2)]; % quats

    % Canonicalize antipodes for SO(3)
    neg = Q(:,1)<0 | (Q(:,1)==0 & (Q(:,2)<0 | (Q(:,2)==0 & (Q(:,3)<0 | (Q(:,3)==0 & Q(:,4)<0)))));
    Q(neg,:) = -Q(neg,:);
    Q = Q ./ vecnorm(Q,2,2); % safety normalize

    t = Q(:,1); x = Q(:,2); y = Q(:,3); z = Q(:,4);
end