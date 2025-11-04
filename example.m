% Some examples illustrating the use of code. First example generates
% B-tensors for powder averaging with a skewed B -tensor. Proceeds via
% isotropic quaternions -> rotations -> B-tensors. Second example is
% standard LTE. Performance is illustrated by computing the variability of
% powder averages over a large set of sample rotations. 
% GFO = geometric filter optimization
% Electro = electrostatic repulsion
%% Example 1: skewed B tensor/SO(3)
B = diag([2 0 3]); % Bmatrix
N = 24; % number of orientations
DDEflag = 0; %no dde symmetry
% cost = 'GFO'; 
cost = 'electro';
Lmax = 8; kappa = 8; s = 8;
Vdeg = util.sobolev(kappa,s,Lmax);

%%% Choose optimization strategy -> output quaternions
% [t, x, y, z, fval, exitflag, output, d] = S3isotropic(N, DDEflag, metric); %local optimizer, electrostatic
[t,x,y,z,fval,exitflag,output] = isotropic_SO3_global(N,cost,Lmax,Vdeg); %global optimizer, cost electrostatic or GFO
% [t,x,y,z] = util.hopf_kronecker(N); %quasi uniform grid

%%% Get Knutsson weights
[a,b,g] = util.quat_to_euler(t,x,y,z); %Euler angles
E = [a,b,g]; %pack
w = knutsson_weights_SO3(E,Lmax,Vdeg); % weights w

%%% Compute rotated B tensors for powder averaging
% First get rotations R from quaternions
Q = cat(2,t,x,y,z);
R = zeros(3,3,N);
for n = 1:N
    R(:,:,n) = util.rotmat(Q(n,:));
end
% Compose rotations to make first rotation the identity - this does not
% change the integral (powder average)
for n = 1:N
    R(:,:,n) = R(:,:,n)*R(:,:,1)';
end
% Brots are rotated B-Tensors
Brots = zeros(3,3,N);
for n = 1:N
    Brots(:,:,n) =  R(:,:,n)*B*R(:,:,n)' ;
end


%% Test powder average
D = diag([0.1 0.1 2]);
nRot = 1e4; %sample rotations
U = util.randrotations(nRot);
S = zeros(nRot,1);
Sw = zeros(nRot,1);
for r = 1:nRot
    Drot = U(:,:,r)*D*U(:,:,r)';
    for n = 1:N
        S(r) = S(r) + exp(-trace(Brots(:,:,n) * Drot));
        Sw(r) = Sw(r) + w(n)*exp(-trace(Brots(:,:,n) * Drot));
    end
    S(r) = S(r)/N;
end
%%% show histograms
tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
ax = nexttile; hold(ax,'on');
ax.FontSize =14;
colors = ax.ColorOrder;
 histogram(S,100,'EdgeColor','none','FaceColor',colors(1,:),'DisplayName','uniform weights','Normalization','pdf');
hold on
histogram(Sw,100,'EdgeColor','none','FaceAlpha',0.7,'FaceColor',colors(2,:),'DisplayName','Knutsson weights',Normalization='pdf');
xlabel('$\bar{S}$',Interpreter='latex'); ylabel('count');
legend;

%% Example 2: GFO for axisymmetric tensors - here LTE

N = 24; % number of orientations
cost = 'GFO';
% cost = 'electro';
Lmax = 8; kappa = 4; s = 5;
Vdeg = util.sobolev(kappa,s,Lmax);

%%%  optimization of directions q
[q, f,flag,out] = S2isotropic(N,cost,Lmax,Vdeg);
q = util.canon_antipodes(q);
% q = util.sph_fibonacci(N,1); %quasi uniform grid

%%% show directions
figure, hold on, sphere(100,'EdgeColor','none'), light, axis equal, view(3)
plot3(q(:,1),q(:,2),q(:,3),'o','MarkerFaceColor','auto');

%%% Get Knutsson weights
w = knutsson_weights_S2(q,Lmax,Vdeg); % weights w
%% Test powder average
b = 1; % b-val
D = diag([0.2 0.2 2]);
nRot = 1e4; %sample rotations
U = util.randrotations(nRot);
S = zeros(nRot,1);
Sw = zeros(nRot,1);
for r = 1:nRot
    Drot = U(:,:,r)*D*U(:,:,r)';
    for n = 1:N
        S(r) = S(r) + exp(-b*q(n,:) * Drot * q(n,:)');
        Sw(r) = Sw(r) + w(n)*exp(-b*q(n,:) * Drot * q(n,:)');
    end
    S(r) = S(r)/N;
end

%%% show histograms
tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
ax = nexttile; hold(ax,'on');
ax.FontSize = 14;
colors = ax.ColorOrder;
histogram(S,100,'EdgeColor','none','FaceColor',colors(1,:),'DisplayName','uniform weights','Normalization','pdf');
hold on
histogram(Sw,100,'EdgeColor','none','FaceAlpha',0.7,'FaceColor',colors(2,:),'DisplayName','Knutsson weights',Normalization='pdf');
xlabel('$\bar{S}$',Interpreter='latex'); ylabel('count');
legend;