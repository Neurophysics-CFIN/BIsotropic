function [q,output]=S2uniform(N,visualize,metric)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses electrostatic repulsion for the sphere in 3 dimensions to find N
% uniformmly distributed directions.
%
% Written by
% Sune Nørhøj Jespersen
% CFIN/MindLab and Dept. of Physics and Astronomy
% Aarhus Universitet
%
% Cell: +45 60896642
% E-mail: sune@cfin.au.dk
% Web: http://www.cfin.au.dk/~sune
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input arguments:
%                   N - integer, number of orientations (pairs)
%                   visuale - logical, to plot results
%                   metric - string, which distance to use, 'geodesic'
%                   (default) or 'euclid'.
%  Output: N directions. Output is information from
%  optimization algorithm

if nargin < 3
    metric = 'geodesic';
end

[x,y,z,fval,exitflag,output,d] = electrostatic_repulsion(N,metric);
% fprintf("Exitflag %s\nObjective function value %f\n",exitflag,fval);
% fprintf("Output structure\n");
% disp(output);

q = cat(2,x,y,z);

if visualize
    scrsz = get(groot,'ScreenSize');
    figure(Position=[scrsz(3)/3 scrsz(4)/3 scrsz(3)/2 scrsz(4)/2]);
    tiledlayout(1,2,'TileSpacing','tight');
    nexttile;
    sphere(100); axis equal; hold on;
    view(-44,18)

    title('q_1')
    for j = 1:N
        x_ = q(j,:);
        plot3(x_(1),x_(2),x_(3),'.r','MarkerSize',12);
        plot3(-x_(1),-x_(2),-x_(3),'.r','MarkerSize',12);
    end
    nexttile;
    axis equal; hold on;
    D = zeros(N,N);
    for i = 1:N
        for j = 1:i-1
            D(i,j) = d([x(i,:),y(i,:),z(i,:)],[x(j,:),y(j,:),z(j,:)]);
        end
    end
    D = .5* (D + D');
    imshow(D,[]);
    colormap parula;
    colorbar
    xlabel('q_i');
    ylabel('q_j');
    title('Normalized distances between directions')
end
end

function [x, y, z, fval, exitflag, output, d] = electrostatic_repulsion(N, metric)
    % Define optimization variables for x, y, z
    x = optimvar('x', N, 1, 'LowerBound', -1, 'UpperBound', 1);
    y = optimvar('y', N, 1, 'LowerBound', -1, 'UpperBound', 1);
    z = optimvar('z', N, 1, 'LowerBound', -1, 'UpperBound', 1);

    % Define the optimization problem
    elecprob = optimproblem;

    % Add the constraint that the points lie on the unit 2-sphere
    elecprob.Constraints.spherec = (x.^2 + y.^2 + z.^2 == 1);

    % Metric selection for distance calculation
    if nargin == 2
        switch metric
            case "euclid"
                d = @(r1, r2) sqrt(max(0,min(1,2-2*abs(r1*r2')))) / 2;
            case "geodesic"
                d = @(r1, r2) acos(max(0,min(1,abs(r1 * r2')))) / pi;
            otherwise
                error("Invalid metric");
        end
    end

    % Define the objective function: energy computation using pairwise distances
    R = [x, y, z]; % Combine the variables into a matrix of coordinates
    energyExpr = fcn2optimexpr(@(X)mean((pdist(X, d) + eps).^-1, 'all'),R);  % Electrostatic energy formula

    % Set the objective function for the optimization problem
    elecprob.Objective = energyExpr;

    % Set up optimization options
    options = optimoptions('fmincon', ...
        'Algorithm', 'interior-point', ...
        'SubproblemAlgorithm', 'cg', ...
        'MaxFunctionEvaluations', 1e3 * N^2, ...
        'MaxIterations', 1e3 * N, ...
        'Display', 'none', ...  % Optional display
        'EnableFeasibilityMode', true);
    x0 = util.sph_fibonacci(N);
    x0 = x0 ./ vecnorm(x0, 2, 2);  % Normalize to lie on the sphere

    init.x = x0(:, 1);
    init.y = x0(:, 2);
    init.z = x0(:, 3);
    % Solve the optimization problem
    [sol, fval, exitflag, output] = solve(elecprob,init, 'Options', options);

    % Extract x, y, z coordinates from the solution
    x = sol.x;
    y = sol.y;
    z = sol.z;
    neg = z < 0 | (z==0 & (y < 0 | (y==0 & (x < 0 )))); %put on one hemisphere
    x(neg) = -x(neg); y(neg) = -y(neg); z(neg) = -z(neg);


    % Display output message if the optimization was not successful
    if exitflag ~= 1
        fprintf("Exitflag %d from %s\n", exitflag, mfilename);
        disp(output.message)
        disp("End")
    end
end
