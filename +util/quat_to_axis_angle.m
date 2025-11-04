function [axis, angle] = quat_to_axis_angle(q)
% Normalize quaternion
q = q / vecnorm(q,2,2);

% Extract quaternion components
x = q(2);
y = q(3);
z = q(4);
t = q(1);
vec = [x; y; z];
% Compute angle of rotation
angle = 2 * atan2(norm(vec),t);

% Compute axis of rotation

axis =  vec/norm(vec);

end
