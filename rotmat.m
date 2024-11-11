function R = rotmat(q)
%By Sune N Jespersen. Email:
% sune@cfin.au.dk Please cite S. Jespersen, MRM XXX (2024) Tested with
% Matlab 2024b.

% Helper function to convert a quaternion to a rotation matrix
%Input: quaternion q = [t,x,y,z]
% Output: Rotation matrix R (3 x 3).

    q = q / norm(q); % Normalize quaternion
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    R = [1 - 2*(q2^2 + q3^2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
         2*(q1*q2 + q0*q3), 1 - 2*(q1^2 + q3^2), 2*(q2*q3 - q0*q1);
         2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1^2 + q2^2)];
end