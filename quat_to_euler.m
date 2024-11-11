function [alpha, beta, gamma] = quat_to_euler(t, x, y, z)
% By Sune N Jespersen. Email: sune@cfin.au.dk
% Please cite S. Jespersen, MRM XXX (2024)
% Tested with Matlab 2024b.
% Convert quaternion (t,x,y,z), with t the scalar part, to Euler angles
% (ZY'Z') of the corresponding rotations.

% Ensure inputs are column vectors
t = t(:);
x = x(:);
y = y(:);
z = z(:);
if any(abs(x.^2 + y.^2 + z.^2 + t.^2 - 1) > 1e-16)
    % warning("Quaternions unnormalized ... normalizing");
    qnorm = sqrt(x.^2 + y.^2 + z.^2 + t.^2);
    % fprintf("Max norm deviation was %g",max(abs(1-qnorm)));
    x = x./qnorm;
    y = y./qnorm;
    z = z./qnorm;
    t = t./qnorm;
end

% Calculate Euler angles
beta = acos(2*(t.^2 + z.^2)-1);
if beta < 1e-16
    beta = 0;
    gamma = 2*atan2(z,t);
    alpha = 0;
elseif abs(beta-pi) < 1e-16
    beta = pi;
    alpha = 0;
    gamma = 2*atan2(-x,y);
else
    thplus = atan2(z,t);
    thminus = atan2(-x,y);
    alpha = (thplus - thminus);
    gamma = (thplus + thminus);
end
end
