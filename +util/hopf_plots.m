function hopf_plots(Q)
% hopf_plots_from_quaternions(Q)
% Q : N x 4 array of quaternions, scalar-first [t x y z].
% Plots:
%   (left) fiber coordinate (S^1) as points on the unit circle
%   (right) base point (S^2) via the Hopf map as points on a sphere
%
% Hopf coordinates from (z1, z2) in C^2:
%   z1 = a + i b, z2 = c + i d  (with q = a + b i + c j + d k)
%   alpha = arg(z1), beta = arg(z2)
%   xi   = alpha + beta          (fiber angle)
%   eta  = alpha - beta
%   theta = 2 * atan2(|z2|, |z1|)
%   n = [sin(theta)cos(eta), sin(theta)sin(eta), cos(theta)] in S^2

    arguments
        Q double {mustBeFinite, mustBeReal, mustBeNonempty, mustBeNcols4(Q)}
    end

    % Normalize quaternions row-wise
    Q = Q ./ vecnorm(Q, 2, 2);

    a = Q(:,1); b = Q(:,2); c = Q(:,3); d = Q(:,4);
    z1 = a + 1i*b;
    z2 = c + 1i*d;

    alpha = atan2(imag(z1), real(z1));
    beta  = atan2(imag(z2), real(z2));

    xi    = mod(alpha + beta, 2*pi);   % fiber angle in [0, 2π)
    eta   = alpha - beta;

    r1    = hypot(a, b);
    r2    = hypot(c, d);
    theta = 2 * atan2(r2, r1);

    % Base point on S^2
    n = [sin(theta).*cos(eta), ...
         sin(theta).*sin(eta), ...
         cos(theta)];

    %---- Plotting ----%
    figure('Color','w'); 

    % Left: fiber S^1
    subplot(1,2,1);
    t = linspace(0, 2*pi, 400);
    plot(cos(t), sin(t), 'k-', 'LineWidth', 0.75); hold on;
    scatter(cos(xi), sin(xi), 18, 'filled');
    axis equal; axis off;
    title('Fiber (S^1)');

    % Right: base S^2
    subplot(1,2,2);
    [Xs, Ys, Zs] = sphere(100);
    surf(Xs, Ys, Zs, 'EdgeColor','none'); hold on;
    plot3(n(:,1), n(:,2), n(:,3), 'o', 'MarkerSize',4, 'MarkerFaceColor',[0 0 0], 'MarkerEdgeColor','none');
    axis equal; view(3);
    camlight headlight; lighting gouraud;
    title('Base (S^2)');

end

% Local validator for 4 columns
function mustBeNcols4(Q)
    if size(Q,2) ~= 4
        error('Input must be an N x 4 array of quaternions [w x y z].');
    end
end
