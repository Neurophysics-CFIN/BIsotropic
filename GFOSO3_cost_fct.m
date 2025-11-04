function val = GFO_cost_fct(x, N, Lmax, Vdeg)
% x: 4N×1 vector → N×4 quaternions (rows). We normalize rows and canonicalize.

Q = reshape(x, [N,4]);
Q = util.normalize_rows(Q);              % enforce unit quaternions
Q = util.canon_antipodes(Q);            % q ≡ -q
nR = size(Q,1);

% Pairwise |dot| using sqrt(s^2) (no abs; stable)
S = Q*Q.';                 % nR×nR
S = sqrt(S.^2);            % entrywise |.| in [0,1]

% Chebyshev U_k(S) recurrence
U_prev = ones(nR);         % U_0
U_curr = 2*S;              % U_1
val = 0;

% Only even ℓ: k = 2ℓ; keep k divisible by 4 → ℓ even
for k = 2 : 2*Lmax
    U_next = 2*S .* U_curr - U_prev;    % U_k
    if mod(k,4) == 0
        ell = k/2;                      % even ℓ
        % Wl  = (1 + (ell*(ell+1))/(kappa^2))^(-s);   % Sobolev weight
        chi = U_next;                                % χ_ℓ = U_{2ℓ}(S)
        sum_all = sum(chi(:));
        term = sum_all / ((2*ell+1)^2);
        val = val + Vdeg(ell/2+1) * term;
    end
    U_prev = U_curr; U_curr = U_next;
end
end