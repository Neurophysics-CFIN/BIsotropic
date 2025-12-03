function val = electroSO3_cost_fct(x, N,addpen)

Q = reshape(x, [N,4]);
%test
% Q = [1 0 0 0;Q];
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
    % val = mean(1 ./ (w(~isnan(w)) + 1e-12)).^2;
end

end

