function energy = electroS2_cost_fct(n, use_proj)
if nargin < 2
    use_proj = false;
end
n = normalize_rows(n);
X   = n*n.';                       % dot products
if use_proj
    X = abs(X);
end
X   = max(-1, min(1, X));
Dch = sqrt(max(0, 2 - 2*X));       % chordal distance on S^2
N   = size(Dch,1);
Dch(1:N+1:end) = Inf;              % remove self
energy = sum( 1./Dch(triu(true(N),1)) );   % unordered pairs once
end


function A = normalize_rows(A)
nrm = sqrt(sum(A.^2, 2));
nrm(nrm == 0) = 1;
A = A ./ nrm;
end