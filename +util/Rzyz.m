function R=Rzyz(a,b,g) %Euler to rot matrix
Rz1=[cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
Ry =[cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
Rz2=[cos(g) -sin(g) 0; sin(g) cos(g) 0; 0 0 1];
R=Rz1*Ry*Rz2;
end