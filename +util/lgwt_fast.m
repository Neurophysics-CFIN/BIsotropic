function [x, w] = lgwt_fast(N, a, b)
%LGWT_FAST  Gauss-Legendre nodes x and weights w on [a,b].
%
%   This is a numerically stable, dependency-free implementation
%   using the Golub–Welsch eigenvalue method.
%
%   INPUT:
%     N : number of nodes (positive integer)
%     a,b : interval endpoints
%
%   OUTPUT:
%     x : N×1 vector of nodes
%     w : N×1 vector of weights
%
%   ACCURACY:
%     ~15 digits for N <= 200.

    if N < 1
        x = []; w = [];
        return;
    end

    % Beta recursion coefficients for Legendre polynomials
    i  = (1:N-1).';
    a0 = zeros(N,1);
    b0 = i ./ sqrt(4*i.^2 - 1);

    % Symmetric tridiagonal Jacobi matrix
    J = diag(b0,1) + diag(a0) + diag(b0,-1);

    % Compute eigenvalues (nodes) and eigenvectors
    [V, D] = eig(J);

    x0 = diag(D);
    [x0, idx] = sort(x0);
    V = V(:, idx);

    % Weights on [-1,1]
    w0 = 2 * (V(1,:).^2).';

    % Map from [-1,1] to [a,b]
    x = (b - a)/2 * x0 + (a + b)/2;
    w = (b - a)/2 * w0;
end
