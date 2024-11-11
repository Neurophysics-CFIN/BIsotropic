function [t, x, y, z, fval, exitflag, output, d] = S3isotropic(N, DDEflag, metric)
% By Sune N Jespersen. Email: sune@cfin.au.dk
% Please cite S. Jespersen, MRM XXX (2024)
% Tested with Matlab 2024b. 

% Input variables:

% - N number of rotations (mandatory)

% - DDEflag (optional), if true (default is false) assumes DDE with rotations of
% q1 = (1,0,0) and q2 = (0,1,0), and takes into account also DDE time
% reversal invariance

% - metric specifies the distance metric, see below. Default is geodesic = angle between two rotations.

% Output variables:

% - t, x, y, z are coordinates of N uniformly distributed directions (unit
% lengh vectors) in 4D. Corresponds to quaterions

% - fval, exitflag and output are optimization diagnostics.
% - d is the metric

% Define optimization variables
t = optimvar('t', N, 1, 'LowerBound', -1, 'UpperBound', 1);
x = optimvar('x', N, 1, 'LowerBound', -1, 'UpperBound', 1);
y = optimvar('y', N, 1, 'LowerBound', -1, 'UpperBound', 1);
z = optimvar('z', N, 1, 'LowerBound', -1, 'UpperBound', 1);


% Define the optimization problem
elecprob = optimproblem;

% Add the constraint that the points lie on the unit sphere in 4D
elecprob.Constraints.spherec = (x.^2 + y.^2 + z.^2 + t.^2 ) == 1;


% Set optimization options
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'EnableFeasibilityMode', true, ...
    'SubproblemAlgorithm', 'cg', ...
    'MaxFunctionEvaluations', 1e3 * N^2, ...
    'MaxIterations', 1e3 * N, ...
    'Display', 'none');  

% Select the metric for distance calculation
if nargin == 3
    switch metric
        case 'euclid'
            d = @(r1, r2) sqrt(sum((r1 - r2 + 1e-16).^2, 2)) / 2;
        case 'geodesic'
            d = @(r1, r2) 2 * acos(sqrt((r1 * r2' + 1e-16).^2))' / pi;
        case 'chordal'
            d = @(r1, r2) 4 * (1- (r1 * r2').^2)' / pi;
        case 'bi-invariant'
            % Axis-angle distance based on rotation angle
            d = @(r1, r2) arrayfun(@(i) axis_angle_distance(r1, r2(i,:)), (1:size(r2, 1))');

        case 'axis-angle'
            % Axis-angle difference, combining rotation angle and axis difference
            d = @(r1, r2) 2 * acos(abs(sum(r1 .* r2, 2))) .* vecnorm(r2 - r1, 2, 2);

        case 'quaternion'
            % Quaternion-based distance, considering axis and angle
            d = @(r1, r2) min(vecnorm(r2 - r1, 2, 2), vecnorm(r2 + r1, 2, 2));

        case 'rotation-axes'
            % Angle between rotation axes combined with rotation angle
            d = @(r1, r2) acos(tensorprod(r1(2:4), r2(:,2:4), 2, 2)' ./ (vecnorm(r1(2:4)) .* vecnorm(r2(:,2:4), 2, 2))) + ...
                2 * acos(abs(sum(r1 .* r2, 2)));

        case 'mod-frobenius'
            % Modified Frobenius norm, difference of rotation matrices
            d = @(r1, r2) arrayfun(@(i) norm(rotmat(r1) - rotmat(r2(i,:)), 'fro'), (1:size(r2, 1))');

        otherwise
            error('Invalid metric');
    end
else
    d = @(r1, r2) 2 * acos(sqrt((r1 * r2' + 1e-16).^2))' / pi; % default metric GEODESIC
end

% Define the objective function

R = [t, x, y, z]; % N x 4
if nargin == 2 && DDEflag
    %adding rotations that leave DDE invariant
    R2 = [R(:,2) - R(:,3), -R(:,1) - R(:,4), R(:,1) - R(:,4), R(:,2) + R(:,3)] / sqrt(2);
    R = cat(1, R, R2); % 2N x 4 matrix
end

% Redefine the line below for other cost functions. Here it is 1/r.
elecprob.Objective = mean(1 ./ (pdist(R, d) + 1e-16));

% Generate an initial guess
x0 = randn(N, 4);
x0 = x0 ./ vecnorm(x0, 2, 2);  % Normalize to lie on the sphere
init.t = x0(:, 1);
init.x = x0(:, 2);
init.y = x0(:, 3);
init.z = x0(:, 4);


% Solve the optimization problem
[sol, fval, exitflag, output] = solve(elecprob, init, 'Options', options);

% Extract the solution
t = sol.t;
x = sol.x;
y = sol.y;
z = sol.z;

end