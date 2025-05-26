function [q1,q2,q3] = quat_to_frame(t, x,y,z)
%given a set of N quaternions (t,x,y,z) with t the scalar part, returns the 
% 3 diffusion wave vectors corresponding to the actions (rotations) of the
% quaternions on the x-, y- and
%z-axes as rows (q are N x 3).
x = reshape(x,[],1);
y = reshape(y,[],1);
z = reshape(z,[],1);
t = reshape(t,[],1);
q1 = [(1 - 2*y.^2 - 2*z.^2), 2*(x.*y + t.*z), 2*(x.*z - t.*y)];
q2 = [2*(x.*y - t.*z), (1 - 2*x.^2 - 2*z.^2), 2*(y.*z + t.*x)];
q3 = cross(q1,q2);
end