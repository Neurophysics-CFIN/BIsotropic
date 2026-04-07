function [x, f,flag,out] = S2isotropic(N, cost, Lmax,bandWeight,x0)
% Settings
% -------------------------------------------------------------------------

%anti-podal flag to identify anti-podal points
use_proj = true;
% fmincon options (interior-point; increase if you want tighter solves)
opt = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','none', ...               % 'iter' if you want logs
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 1e5, ...
    'StepTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6);

% Initial step: use spherical Fibonacci as a good seed
if ~(nargin >= 5 && ~isempty(x0))
    x0 = util.sph_fibonacci(N);
end
% x0 = x0 + randn(size(x0))*0.08;
% x0 = util.normalize_rows(x0);

if strcmp(cost,'electro') % Electrostatic
    obj = @(x)electroS2_cost_fct(x,use_proj);
elseif strcmp(cost,'GFO')
    obj = @(x)GFOS2_cost_fct(x,Lmax,bandWeight);
end

lb = -ones(size(x0));
ub = ones(size(x0));
problem = struct('objective', obj, ...
    'x0', x0, ...
    'Aineq', [], 'bineq', [], ...
    'Aeq', [], 'beq', [], ...
    'lb', lb, 'ub', ub, ...
    'nonlcon', [], ...
    'solver', 'fmincon', ...
    'options', opt);
[x, f,flag,out] = fmincon(problem);
x = util.normalize_rows(x);
x = util.canon_antipodes(x);
end


