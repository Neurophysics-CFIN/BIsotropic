function [w, meta] = knutsson_weights_SO3(E,Lmax,bandWeight,meta)

if nargin == 3
    % build_B_complex: returns complex B and the degree l for each row
    [Bc, N, ~, l_of_row] = build_B(E, Lmax);   % complex
    W = bandWeight(l_of_row/2+1)';
    W    = [W; W];     
    % Build real-augmented system:
    A   = [real(Bc); imag(Bc)];                        % (2M x N)
    g0c = zeros(size(Bc,1),1); g0c(1) = 1;             % complex target (only l=0,m=n=0)
    b   = [real(g0c); imag(g0c)];                      % real target (2M x 1)
    % duplicate for Re/Im
    M = size(Bc,1);
    W(M+1) = 0;
elseif  nargin ==4 & ~isempty(meta)
    N=size(E,1);
    Bc = meta.Bc;
    l_of_row = meta.l_of_row;
    W = meta.W;
    A = meta.A;
    b = meta.b;
else
    error("Syntax error");
end
% lambda = 1e-6;                                     % Tikhonov regularization
AtA = A.'*(W.*A);% + lambda*eye(N);
lam = 1e-6 * trace(AtA) / max(N,1); %ridge
AtA =  AtA + lam*eye(N);
rhs = A.'*(W.*b);
w   = AtA\rhs;
w   = w/sum(w);     % normalize to sum=1
if nargout == 2
    meta.Bc = Bc;
    meta.l_of_row = l_of_row;
    meta.V = W;           % weights used in solve (2M x 2M)
    meta.A = A;
    meta.b = b;
end
end

function [B, N, M, l_of_row] = build_B(E,Lmax)
N = size(E,1);
ells = 0:2:Lmax;  % even l only due to symmetry (180 around eigenvectors)
% rows = sum( (2*ells + 1).^2 );
rows = sum( 2*ells.*(ells+1) + 1 );
B = zeros(rows, N);
l_of_row = zeros(rows,1);
r = 1;
for l = 0:2:Lmax
    for m = 0:l
        for n = -l:l
            if (m > 0) || (m == 0 && n >= 0) %accounting for signal reality by only including independent coefficients
                for i = 1:N
                    a = E(i,1); b = E(i,2); g = E(i,3);
                    B(r,i) = conj(util.wignerD(l,m,n,a,b,g));
                end
                l_of_row(r) = l;
                r = r + 1;
            end
        end
    end
end
M = rows;
end