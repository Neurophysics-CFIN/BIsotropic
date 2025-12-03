function q = rotm2quat2(R)
    % ROTM2QUAT Converts a rotation matrix to a quaternion
    % Input: R - 3x3 rotation matrix
    % Output: q - 1x4 quaternion [t, x, y, z]

    % Ensure the matrix is a valid rotation matrix
    if ~isequal(size(R), [3 3])
        error('Input matrix must be 3x3.');
    end
    
    % Check for the orthonormality of the matrix
    if norm(R*R' - eye(3), 'fro') > 1e-6 || abs(det(R) - 1) > 1e-6
        error('Input matrix must be a valid rotation matrix.');
    end
    
    % Calculate the trace of the matrix
    t = trace(R);
    
    % Initialize quaternion
    q = zeros(1, 4);
    
    if t > 0
        % Positive trace case
        r = sqrt(1 + t);
        q(1) = 0.5 * r;  % t
        q(2) = (R(3,2) - R(2,3)) / (2 * r);  % x
        q(3) = (R(1,3) - R(3,1)) / (2 * r);  % y
        q(4) = (R(2,1) - R(1,2)) / (2 * r);  % z
    else
        % Handle the edge cases based on which diagonal term is the largest
        if (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
            % R(1,1) is the largest
            r = sqrt(1 + R(1,1) - R(2,2) - R(3,3));
            q(2) = 0.5 * r;  % x
            q(1) = (R(3,2) - R(2,3)) / (2 * r);  % t
            q(3) = (R(1,2) + R(2,1)) / (2 * r);  % y
            q(4) = (R(1,3) + R(3,1)) / (2 * r);  % z
        elseif R(2,2) > R(3,3)
            % R(2,2) is the largest
            r = sqrt(1 + R(2,2) - R(1,1) - R(3,3));
            q(3) = 0.5 * r;  % y
            q(1) = (R(1,3) - R(3,1)) / (2 * r);  % t
            q(2) = (R(1,2) + R(2,1)) / (2 * r);  % x
            q(4) = (R(2,3) + R(3,2)) / (2 * r);  % z
        else
            % R(3,3) is the largest
            r = sqrt(1 + R(3,3) - R(1,1) - R(2,2));
            q(4) = 0.5 * r;  % z
            q(1) = (R(2,1) - R(1,2)) / (2 * r);  % t
            q(2) = (R(1,3) + R(3,1)) / (2 * r);  % x
            q(3) = (R(2,3) + R(3,2)) / (2 * r);  % y
        end
    end
    
    % Normalize the quaternion to ensure it is valid
    q = q / norm(q);
end
