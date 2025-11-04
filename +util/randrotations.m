function U = randrotations(K)
% Haar random rotations via random quaternions
U = zeros(3,3,K);
q = randn(K,4); q = q./vecnorm(q,2,2);
for k = 1:K
    U(:,:,k) = util.quat_to_R(q(k,:));
end
end