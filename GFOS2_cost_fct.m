function Fsum  = GFOS2_cost_fct(n, Lmax, bandWeight)
n    = util.normalize_rows(n);
X    = n*n.';                           % x_ij = n_i · n_j
ells = 0:2:Lmax;
bandWeight = bandWeight(:).'; bandWeight(1) = 0;
Fsum = 0;
for li = 2:numel(ells)                  % skip ℓ=0
    l     = ells(li);
    Pl    = legendreP_matrix(l, X);     % NxN, element-wise P_l(x)
    Kfac  = (2*l+1)/(4*pi);             % addition-theorem factor
    Fsum  = Fsum + bandWeight(li) * Kfac * sum(Pl, 'all');
end

end

function X = legendreP_matrix(l, Xdot)
% Return P_l(Xdot), element-wise (Xdot is NxN of dot products)
% Uses MATLAB's LEGENDRE for scalar x via vectorization over unique values
x = Xdot(:);
% Handle endpoints robustly
x(x> 1) = 1;
x(x<-1) = -1;
% P = legendre(l, x);         % returns (m+1) x numel(x); take m=0 row
% P0 = P(1,:).';
P0 = util.legendreP(l, x)';         % returns 1 x numel(x); 
X = reshape(P0, size(Xdot));
end
