function n = sph_fibonacci(N, use_proj)
% Spherical (or projective-hemisphere) Fibonacci
% If use_proj==true, generate N points on the UPPER hemisphere (z>0)
% so ±n never both appear at init. Otherwise, full sphere.

if nargin < 2, use_proj = false; end

i   = (0:N-1)';
phi = (1 + sqrt(5))/2;          % golden ratio
theta = 2*pi*i/phi;             % azimuths

if ~use_proj
    % Full sphere: z in (-1,1)
    z = -1 + 2*(i + 0.5)/N;
else
    % Projective-safe: UPPER hemisphere only (z in (0,1))
    % This avoids antipodal duplicates at init.
    z = (i + 0.5)/N;            % no 0 or 1 exactly; no equator ties
end

r = sqrt(max(0, 1 - z.^2));
n = [r.*cos(theta), r.*sin(theta), z];
n = util.normalize_rows(n);

% Optional nano-jitter to avoid any accidental coincidences for tiny N
if use_proj
    n = util.normalize_rows(n + 1e-12*randn(size(n)));
end
end