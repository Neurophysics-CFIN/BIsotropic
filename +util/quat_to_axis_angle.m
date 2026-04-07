function [axis, angle] = quat_to_axis_angle(q)
%QUAT_TO_AXIS_ANGLE  Convert a quaternion to axis–angle representation.
%
%   [axis, angle] = QUAT_TO_AXIS_ANGLE(q)
%
%   Converts a rotation represented by a quaternion into its equivalent
%   axis–angle form. The quaternion is normalized internally before the
%   conversion.
%
%   INPUT
%     q : 4-by-1 quaternion in the convention [t x y z], where t is the
%         scalar part. The quaternion need not be unit length.
%
%   OUTPUTS
%     axis  : 3-by-1 unit vector giving the rotation axis.
%     angle : Rotation angle (in radians).
%
%   DETAILS
%     - The quaternion is first normalized to unit length.
%     - The rotation angle is computed as
%           angle = 2 * atan2(norm([x y z]), t).
%     - The rotation axis is obtained by normalizing the vector part
%       [x y z]. This assumes a nonzero rotation angle.
%
%   NOTES
%     - For quaternions representing a zero rotation (vector part equal to
%       zero), the rotation axis is undefined; this case is not explicitly
%       handled in the code.
%     - The returned angle lies in [0, 2*pi).
%
%   EXAMPLE
%     q = [cos(pi/4); sin(pi/4); 0; 0];   % 90-degree rotation about x-axis
%     [axis, angle] = quat_to_axis_angle(q);

% Normalize quaternion
q = q / vecnorm(q,2,2);

% Extract quaternion components
t = q(1);
x = q(2);
y = q(3);
z = q(4);

vec = [x; y; z];
% Compute angle of rotation
angle = 2 * atan2(norm(vec),t);

% Compute axis of rotation

axis =  vec/norm(vec);

end
