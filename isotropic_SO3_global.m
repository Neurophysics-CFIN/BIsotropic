function [t,x,y,z,fval,exitflag,output] = isotropic_SO3_global(N, cost, Lmax, bandWeight)
% Global SO(3) design 
% Global stage: particleswarm (bounds only)  → HybridFcn: fmincon polish
% Input:
%   N = number of rotations; cost = cost function, "electro" or "GFO"; Lmax = maximum ℓ for GFO; Vdeg = degree
%   weighting as function of ℓ (ℓ/2 + 1 vector)
% Returns unit quaternions [t x y z] (Nx1 each). t is scalar part
%
% Usage example:
% [t,x,y,z,fval,exitflag,output] = isotropic_SO3_global(24,"GFO",4,[1,1/2,1/8])


dim = 4*N;                 % optimization variables are stacked [t;x;y;z] (4N×1)
lb  = -ones(dim,1);        % PSO bounds (loose; we renormalize inside the objective)
ub  =  ones(dim,1);

% --- Deterministic initializer (Hopf–Kronecker on S³, antipodal-canonical) ---
[t,x,y,z] = util.hopf_kronecker(N);
x0 = cat(1,t(:),x(:),y(:),z(:));

if strcmp(cost,"GFO")
    % --- Objective handle (numeric, bounds-only; normalizes inside) ---
    obj = @(x) GFOSO3_cost_fct(x,N,Lmax,bandWeight);
elseif strcmp(cost,"electro")
    obj = @(x) electroSO3_cost_fct(x,N,false);
else
    error('Invalid cost function specified. Choose either "GFO" or "electro".');
end

% --- particleswarm options (seed swarm near Hopf–Kronecker, enable hybrid polish) ---
% Build an initial swarm around x0 with small random jitters
swarmSz = max(40, 20*N);                         % heuristic
J = 0.05;                                        % jitter amplitude (in [-J,J])
InitSwarm = repmat(x0.', swarmSz, 1) + J*(2*rand(swarmSz,dim)-1);

% HybridFcn: fmincon polish on same bounds; no nonlinear constraints needed
hybOpt = optimoptions('fmincon', 'Algorithm','interior-point', 'Display','none', ...
    'SpecifyObjectiveGradient',false,'MaxFunctionEvaluations',1000);

psopt = optimoptions('particleswarm', ...
    'SwarmSize', swarmSz, ...
    'MaxIterations', 1000, ...
    'MaxStallIterations', 60, ...
    'UseParallel', false, ...            % set true if you have Parallel Toolbox
    'Display', 'off', ...
    'InitialPoints', InitSwarm, ...
    'HybridFcn', {@fmincon, hybOpt});


% --- Run particleswarm ---
[xbest, ~, ~, outPS] = particleswarm(obj, dim, lb, ub, psopt);

% --- Optional final polish (again) with tighter tolerances ---
if strcmp(cost,"electro")
    obj = @(x) electroSO3_cost_fct(x,N,false);
end
polishOpt = optimoptions('fmincon','Algorithm','interior-point','SubproblemAlgorithm',...
    'cg', 'Display','none','MaxIterations', 5000, 'OptimalityTolerance', 1e-6,...
    'StepTolerance',1e-10, 'MaxFunctionEvaluations',1e5);
problem = createOptimProblem('fmincon','objective',obj,'x0',xbest,'lb',lb,'ub',ub,'options',polishOpt);
[xbest, fbest, eflagPolish, outPolish] = fmincon(problem);

% --- Extract unit quaternions from best vector (normalize + canonicalize) ---
Q = reshape(xbest, [N,4]);
Q = util.normalize_rows(Q);
Q = util.canon_antipodes(Q);

t = Q(:,1); x = Q(:,2); y = Q(:,3); z = Q(:,4);
fval = fbest;
exitflag = eflagPolish;
output = struct('particleswarm', outPS, 'polish', outPolish);

end
