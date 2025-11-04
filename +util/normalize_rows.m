function A = normalize_rows(A)
nrm = sqrt(sum(A.^2,2)); nrm(nrm==0) = 1; A = A./nrm;
end
