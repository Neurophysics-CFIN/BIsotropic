function [E, w] = euler_so3_grid(Lmax)
% Exact band-limited SO(3) quadrature up to bandwidth Lmax.
% Produces a direct-product Euler grid (alpha,beta,gamma)
% with weights w(j) so that integrals up to l<=Lmax are exact.
% N = (2Lmax +1)^2 * (Lmax +1)
% Output:
%   E : [N × 3] array of [alpha, beta, gamma]
%   w : [1 × N] weights, sum(w) = 1

    % ---- 1) alpha, gamma: uniform Fourier grids ----
    Na = 2*Lmax + 1;              % uniform samples for exact Fourier modes
    Ng = 2*Lmax + 1;
    alpha = (0:Na-1) * 2*pi/Na;   % row
    gamma = (0:Ng-1) * 2*pi/Ng;

    % ---- 2) beta: Gauss–Legendre quadrature for cos(beta) ----
    % Integrates polynomials up to degree 2Lmax exactly.
    Nbeta = Lmax + 1;
    [x, w_leg] = util.lgwt_fast(Nbeta, -1, 1);     % x = cos(beta)
    beta = acos(x);                      % convert to beta angles

    % ---- 3) Tensor product grid ----
    [A, B, G] = ndgrid(alpha, beta, gamma);

    % ---- 4) Weights: Haar measure = (1/8π²) dα sinβ dβ dγ ----
    % Combined product weights:
    %
    %  w = w_α * w_γ * (w_leg * sin(beta))
    %
    % Normalize so sum(w)=1 
    %
    w_alpha  = (1/Na) * ones(1,Na);
    w_gamma  = (1/Ng) * ones(1,Ng);
    w_beta   = w_leg(:)' ./ 2;      % because ∫_{-1}^1 dx = ∫_0^π sinβ dβ

    % Expand weights to full 3D grid  
    W = bsxfun(@times, ...
            reshape(w_alpha, [Na 1 1]), ...
            bsxfun(@times, reshape(w_beta, [1 Nbeta 1]), ...
                          reshape(w_gamma, [1 1 Ng])));

    % Normalize to sum(w)=1
    W = W / sum(W(:));

    % ---- 5) Flatten ----
    E = [A(:), B(:), G(:)];
    w = W(:).';
end
