function quaternions = S2thenSO3optim(N,metric)
By Sune N Jespersen. Email:
% sune@cfin.au.dk Please cite S. Jespersen, MRM XXX (2024) Tested with
% Matlab 2024b.

% 2-step algorithm for uniform rotations (as quaternions, see also
% 'S3isotropic.m'), exploiting the decomposition of SO(3) into S2 x S1.
% First the axis of rotation is selected by 'ordinary electrostatic
% repulsion' on the sphere (see Jones et al, MRM 1999;42(3), p515). Second,
% with fixed rotation axes,  the rotation angles are determined by
% minimizing the overall 'energy' on SO(3). 

% See 'S3isotropic.m' for input and output arguments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Generate N directions on the unit sphere
% directions = optimize_directions_on_sphere(N);
directions = S2uniform(N,false);
% Step 2: Minimize electrostatic energy of quaternions while keeping directions fixed
if nargin == 1
    metric = 'geodesic';
end
quaternions = optimize_quaternions(directions,metric);

end

function directions = optimize_directions_on_sphere(N)
% Initial guess: random points on the sphere
directions = randn(N, 3);
directions = directions ./ vecnorm(directions, 2, 2);  % Normalize to unit length

% Apply constrained optimization to keep points on the sphere
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');

directions = fmincon(@(x) coulomb_energy(x, N), directions(:), [], [], [], [], [], [], @(x) sphere_constraint(x, N), options);
directions = reshape(directions, [N, 3]);
directions = directions ./ vecnorm(directions,2,2);
end

function energy = coulomb_energy(directions, N)
directions = reshape(directions, [N, 3]);
S2_dist = @(r1, r2) acos(max(-1,min(1,r1 * r2'))) / pi; %Geodesic
% S2_dist = @(r1, r2) 1./sqrt(vecnorm(r1-r2,2,2));  %Eucledian
energy = mean((pdist(directions, S2_dist) + eps).^-1);
% energy = mean(log(pdist(directions, S2_dist) + eps)); %Real Coulomb in 2D
end

function [c, ceq] = sphere_constraint(directions, N)
% Constraint to keep points on the 3D sphere
directions = reshape(directions, [N, 3]);
ceq = sum(directions.^2, 2) - 1;  % Ensure each direction has norm 1
c = [];
end

function q = direction_to_quaternion(direction, angle)
% Given a direction, compute the quaternion that rotates the z-axis onto it
% and then rotates around that direction by an angle (second degree of freedom)
q = zeros(length(angle),4);
z_axis = [0, 0, 1];
theta = zeros(length(angle),1);
N = length(angle);
% Step 1: Rotate z-axis to the direction
axis = cross(repmat(z_axis, N, 1), direction);
axis_norm = vecnorm(axis,2,2);
ix = axis_norm == 0;
% If the direction is along z, just use the angle rotation around z-axis
q(ix,:) = [cos(angle(ix) / 2), zeros(length(find(ix)),2), sin(angle(ix) / 2)];
axis(~ix,:) = axis(~ix,:) ./ axis_norm(~ix,:);  % Normalize the axis
theta(~ix) = acos(dot(repmat(z_axis, sum(~ix), 1), direction(~ix, :), 2));

q1 = [cos(theta / 2), axis .* sin(theta / 2)];  % Quaternion to align z-axis with direction

% Step 2: Rotate around the new axis by the given angle
q2 = [cos(angle / 2), direction .* sin(angle / 2)];  % Quaternion for rotation around direction

% Combine the two quaternions
q(~ix,:) = quat_multiply(q2(~ix,:), q1(~ix,:));
end

function q = quat_multiply(q1, q2)
    % Quaternion multiplication for two arrays of quaternions
    % q1 and q2 are both N x 4 arrays, where each row is a quaternion [w, x, y, z]
    % Returns an N x 4 array of the products of q1 and q2.
    
    % Extract the scalar and vector parts of each quaternion
    w1 = q1(:, 1); x1 = q1(:, 2); y1 = q1(:, 3); z1 = q1(:, 4);
    w2 = q2(:, 1); x2 = q2(:, 2); y2 = q2(:, 3); z2 = q2(:, 4);
    
    % Compute the product for all quaternions using vectorized operations
    w = w1 .* w2 - x1 .* x2 - y1 .* y2 - z1 .* z2;
    x = w1 .* x2 + x1 .* w2 + y1 .* z2 - z1 .* y2;
    y = w1 .* y2 - x1 .* z2 + y1 .* w2 + z1 .* x2;
    z = w1 .* z2 + x1 .* y2 - y1 .* x2 + z1 .* w2;
    
    % Combine the scalar and vector parts back into the quaternion form
    q = [w, x, y, z];
end

function quaternions = optimize_quaternions(directions,metric)
% Optimize quaternion configuration while keeping directions fixed
% options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');
N = size(directions,1);
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'SubproblemAlgorithm', 'cg', ...
    'MaxFunctionEvaluations', 1e3 * N^2, ...
    'MaxIterations', 1e3 * N, ...
    'Display', 'none', ...  % Optional display
    'EnableFeasibilityMode', true);
% Optimize only the angles around each direction
angles = rand(size(directions, 1), 1) * 2 * pi;  % Initial random guess for angles
lb = zeros(size(angles)); % Lower bound for angles (0)
ub = 2 * pi * ones(size(angles)); % Upper bound for angles (2*pi)
angles = fmincon(@(angles) quaternion_energy(angles, directions, metric), angles, [], [], [], [], lb, ub,[], options);

% Recompute quaternions with the optimized angles
quaternions = direction_to_quaternion(directions, angles);

% Normalize quaternions to lie on the 4D unit sphere
quaternions = quaternions ./ vecnorm(quaternions,2,2);
end

function energy = quaternion_energy(angles, directions,metric)
if nargin ==3 
    switch metric
            case 'euclid'
                q_dist = @(r1, r2) sqrt(min(sum((r1 - r2 + 1e-16).^2, 2),sum((r1 + r2 + 1e-16).^2, 2))) / 2;
            case 'geodesic'
                q_dist = @(r1, r2) 2 * acos(abs(r1 * r2'))' / pi;
            case 'chordal'
                q_dist = @(r1, r2) 4 * (1- (r1 * r2').^2)' / pi;
            case 'bi-invariant'
                % Axis-angle distance based on rotation angle
                q_dist = @(r1, r2) arrayfun(@(i) axis_angle_distance(r1, r2(i,:)), (1:size(r2, 1))');
                
            case 'axis-angle'
                % Axis-angle difference, combining rotation angle and axis difference
                q_dist = @(r1, r2) 2 * acos(abs(sum(r1 .* r2, 2))) .* vecnorm(r2 - r1, 2, 2);
                
            case 'quaternion'
                % Quaternion-based distance, considering axis and angle
                q_dist = @(r1, r2) min(vecnorm(r2 - r1, 2, 2), vecnorm(r2 + r1, 2, 2));
                
            case 'rotation-axes'
                % Angle between rotation axes combined with rotation angle
                q_dist = @(r1, r2) acos(tensorprod(r1(2:4), r2(:,2:4), 2, 2)' ./ (vecnorm(r1(2:4)) .* vecnorm(r2(:,2:4), 2, 2))) + ...
                              2 * acos(abs(sum(r1 .* r2, 2)));
                
            case 'mod-frobenius'
                % Modified Frobenius norm, difference of rotation matrices
                q_dist = @(r1, r2) arrayfun(@(i) norm(rotmat(r1) - rotmat(r2(i,:)), 'fro'), (1:size(r2, 1))');
                
            otherwise
                error('Invalid metric');
    end
else
    q_dist = @(r1, r2) 2 * acos(abs(r1 * r2'))' / pi;
end
% q_dist = @(r1, r2) sqrt(sum((r1 - r2 + 1e-16).^2, 2)) / 2;
% Calculate electrostatic energy of quaternions with varying angles

q_mat1 = direction_to_quaternion(directions, angles);
q_mat2 = [q_mat1(:,2) - q_mat1(:,3), -q_mat1(:,1) - q_mat1(:,4), q_mat1(:,1) - q_mat1(:,4), q_mat1(:,2) + q_mat1(:,3)] / sqrt(2);
q_mat = cat(1,q_mat1,q_mat2);
[q1,q2] = quat_to_2dframe(q_mat(:,1),q_mat(:,2),q_mat(:,3),q_mat(:,4));
Q = cat(2,q1,q2);
q_dist = @(r1, r2) abs(tensorprod(r1(1:3),r2(:,4:6),2,2)) + abs(tensorprod(r1(4:6),r2(:,1:3),2,2) + tensorprod(r1(1:3),r2(:,1:3),2,2) + tensorprod(r1(4:6),r2(:,4:6),2,2)  );
energy = mean(1 ./ (pdist(Q, q_dist) + eps)); % Add small value to avoid singularities
% energy = mean(1 ./ (pdist(q_mat, q_dist) + eps)); % Add small value to avoid singularities
end

% Helper function to compute axis-angle distance between two quaternions
function dist = axis_angle_distance(q1, q2)
    % Convert the quaternions to rotation matrices
    R1 = rotmat(q1);
    R2 = rotmat(q2);
    
    % Compute the relative rotation matrix
    R_rel = R1' * R2;
    
    % Compute the trace of the relative rotation matrix
    trace_R = trace(R_rel);
    trace_R = max(min(trace_R, 3), -1);%
    
    % Compute the rotation angle (using the trace formula)
    theta = acos((trace_R - 1) / 2);
    
    % Return the rotation angle as the distance
    dist = theta;
end