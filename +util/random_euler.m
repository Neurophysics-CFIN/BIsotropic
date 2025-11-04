function E = random_euler(N)
alpha=2*pi*rand(N,1); gamma=2*pi*rand(N,1); v=2*rand(N,1)-1; beta=acos(v);
E=[alpha,beta,gamma];
end
