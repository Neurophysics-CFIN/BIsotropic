function R = quat_to_R(q)
%QUAT_TO_R Convert quaternion(s) to rotation matrix/matrices.
%   q: [N x 4] or [1 x 4], scalar-first [w x y z]
%   R: [3 x 3 x N]

    if size(q,2) ~= 4
        error('q must be N×4 with [w x y z].');
    end

    % If input is 1x4, this still works with N = 1.
    w = q(:,1);
    x = q(:,2);
    y = q(:,3);
    z = q(:,4);

    % Optional: normalize quaternions to be safe
    % s = sqrt(w.^2 + x.^2 + y.^2 + z.^2);
    % w = w./s; x = x./s; y = y./s; z = z./s;

    N = size(q,1);
    R = zeros(3,3,N,class(q));  % preserve numeric type

    R(1,1,:) = 1 - 2*(y.^2 + z.^2);
    R(1,2,:) = 2*(x.*y - z.*w);
    R(1,3,:) = 2*(x.*z + y.*w);

    R(2,1,:) = 2*(x.*y + z.*w);
    R(2,2,:) = 1 - 2*(x.^2 + z.^2);
    R(2,3,:) = 2*(y.*z - x.*w);

    R(3,1,:) = 2*(x.*z - y.*w);
    R(3,2,:) = 2*(y.*z + x.*w);
    R(3,3,:) = 1 - 2*(x.^2 + y.^2);
end
