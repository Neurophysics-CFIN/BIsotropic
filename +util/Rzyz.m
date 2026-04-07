function R = Rzyz(a,b,g) 
%Euler angles (Rz(a)Ry(b)Rz(g)) to rotation matrix
a = a(:);
b = b(:);
g = g(:);
R = zeros(3,3,length(a));
for j = 1:length(a)
    Rz1 = [cos(a(j)) -sin(a(j)) 0; sin(a(j)) cos(a(j)) 0; 0 0 1];
    Ry = [cos(b(j)) 0 sin(b(j)); 0 1 0; -sin(b(j)) 0 cos(b(j))];
    Rz2 = [cos(g(j)) -sin(g(j)) 0; sin(g(j)) cos(g(j)) 0; 0 0 1];
    R(:,:,j) = Rz1*Ry*Rz2;
end
end