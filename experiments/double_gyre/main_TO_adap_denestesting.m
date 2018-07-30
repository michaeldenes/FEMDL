addpath('../../src'); clear all; colormap jet; clc

%% flow map
t0 = 0; tf = 1; nt = 6; tspan = linspace(t0,tf,nt); 
T = @(x) flow_map(@double_gyre,x,tspan);

%% data points
n = 25; x = linspace(0,1,n);                            
[X,Y] = meshgrid(x,x); p{1} = [X(:) Y(:)];
% p{1} = rand(n^2,2); 
% load p_Meiss
% p1 = rand(20000,2)*diag([1, 0.5])+ones(20000,2)*diag([0 0.5]);
% p2 = rand(200,2)*diag([1, 0.5]);
% p{1} = [p1; p2]; n = size(p{1},1);
pb = [1:n^2; 1:n^2]';


%% time integration
tic; P = T(p{1}); toc
for k = 1:nt, p{k} = [P(:,k) P(:,k+nt)]; end;
%% time integration

p{1}(2,:)
p{2}(2,:)
p{3}(2,:)
p{4}(2,:)
p{5}(2,:)
p{6}(2,:)
%% assembly (for missing data case)
tic; pm = 1;                                % percentage of nodes to remove
% pm = 0.4;                                   % for missing data case
D = sparse(n^2,n^2); M = sparse(n^2,n^2); 
for k = 1:nt
    r = randperm(n^2,floor(pm*n^2))'; 
    tr = delaunay(p{k}(r,:)); 
    t = [r(tr(:,1)), r(tr(:,2)), r(tr(:,3))];
    A = kron([1 0 1],ones(size(t,1),1));      % 2 x 2 identity matrix
    [Dt{k},Mt{k}] = assemble(p{k},t,pb,A); 
    D = D + Dt{k};  M = M + Mt{k}; 
end;
toc


%% some plots
figure(1); scatter(p{1}(100:125,1), p{1}(100:125,2),15, p{1}(100:125,1));
xlabel('$x$'); ylabel('$y$');

%% some plots
figure(2); scatter(p{2}(:,1), p{2}(:,2),3, p{1}(:,1));
xlabel('$x$'); ylabel('$y$');


%% some plots
figure(3); scatter(p{3}(:,1), p{3}(:,2),3, p{1}(:,1));
xlabel('$x$'); ylabel('$y$');



%% triangulation some plots
t = delaunay(p{1}); 
figure(1); triplot(t, p{1}(:,1), p{1}(:,2));
xlabel('$x$'); ylabel('$y$');

%% triangulation plots
t = delaunay(p{2}); 
figure(2); triplot(t, p{2}(:,1), p{2}(:,2));
xlabel('$x$'); ylabel('$y$');


%% triangulation some plots
t = delaunay(p{3}); 
figure(3); triplot(t, p{3}(:,1), p{3}(:,2));
xlabel('$x$'); ylabel('$y$');

%% Testing assemble.m
%% ASSEMBLE stiffness and mass matrices
%
% [D,M] = ASSEMBLE(p,t,pb,G) computes the stiffness matrix D and the mass matrix M
%   p: (n x 2), one node per row
%   t: (m x 3), integers, each row defines a triangle by indexing into p
%   pb: (n x 2), node pb(i,2) maps to pb(i,1) (for perodic boundaries)
%   G: (m x 3), each row defines a tensor on the correspondig triangle
%
% based on code from ifem by Long Chen
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 
G = A
n = max(pb(:,2)); m = size(t,1);

[dphi,area] = gradbasis(p{1},delaunay(p{1}));

% assembly
D = sparse(n,n); M = sparse(n,n);
for i = 2:2
    for j = 2:2

        Dij = -area.*(dphi(:,1,i).*G(:,1).*dphi(:,1,j) ...
                    + dphi(:,1,i).*G(:,2).*dphi(:,2,j) ...
                    + dphi(:,2,i).*G(:,2).*dphi(:,1,j) ...
                    + dphi(:,2,i).*G(:,3).*dphi(:,2,j));
        Mij = area/12.*ones(size(dphi,1),1);
        I = pb(t(:,i),2); J = pb(t(:,j),2);
        if (j==i)
            D = D + sparse(I,J,Dij,n,n);
            M = M + sparse(I,J,Mij+area/12,n,n);
        else
            D = D + sparse([I;J],[J;I],[Dij; Dij],n,n);   
            M = M + sparse([I;J],[J;I],[Mij; Mij],n,n);
        end        
    end
end
%% Testing gradbasis.m
tr = delaunay(p{1});

v(:,:,1) = p{1}(tr(:,3),:)-p{1}(tr(:,2),:);
v(:,:,2) = p{1}(tr(:,1),:)-p{1}(tr(:,3),:);
v(:,:,3) = p{1}(tr(:,2),:)-p{1}(tr(:,1),:);
area = 0.5*(-v(:,1,3).*v(:,2,2) + v(:,2,3).*v(:,1,2));
dphi(:,:,:) = [-v(:,2,:)./(2*area), v(:,1,:)./(2*area)];
tester = -v(:,2,:)./(2*area);
area = abs(area);



%% ASSEMBLE stiffness and mass matrices
%
% [D,M] = ASSEMBLE(p,t,pb,G) computes the stiffness matrix D and the mass matrix M
%   p: (n x 2), one node per row
%   t: (m x 3), integers, each row defines a triangle by indexing into p
%   pb: (n x 2), node pb(i,2) maps to pb(i,1) (for perodic boundaries)
%   G: (m x 3), each row defines a tensor on the correspondig triangle
%
% based on code from ifem by Long Chen
%
% (C) 2017 by O. Junge and G. Froyland, see COPYRIGHT 
G = A;
n = max(pb(:,2)); m = size(t,1);

[dphi,area] = gradbasis(p{1},delaunay(p{1}));

% assembly
D = sparse(n,n); M = sparse(n,n);
for i = 1:1
    for j = i:1
        Dij = -area.*(dphi(:,1,i).*G(:,1).*dphi(:,1,j) ...
                    + dphi(:,1,i).*G(:,2).*dphi(:,2,j) ...
                    + dphi(:,2,i).*G(:,2).*dphi(:,1,j) ...
                    + dphi(:,2,i).*G(:,3).*dphi(:,2,j));
        Mij = area/12.*ones(size(dphi,1),1);
        I = pb(t(:,i),2); J = pb(t(:,j),2);
        if (j==i)
            D = D + sparse(I,J,Dij,n,n);
            M = M + sparse(I,J,Mij+area/12,n,n)
        else
            D = D + sparse([I;J],[J;I],[Dij; Dij],n,n);   
            M = M + sparse([I;J],[J;I],[Mij; Mij],n,n);
        end        
    end
end


%% Testing gradbasis
p0 = p{1}
t0 = delaunayTriangulation(p{1})

v(:,:,1) = p0(t0(:,3),:)-p0(t0(:,2),:);
v(:,:,2) = p0(t0(:,1),:)-p0(t0(:,3),:);
v(:,:,3) = p0(t0(:,2),:)-p0(t0(:,1),:);
area = 0.5*(-v(:,1,3).*v(:,2,2) + v(:,2,3).*v(:,1,2));
dphi(:,:,:) = [-v(:,2,:)./(2*area), v(:,1,:)./(2*area)];
area = abs(area);



%% remove all zero rows and columns
S = sum(abs(D));
I = find(abs(S)>eps);
D = D(I,I); M = M(I,I); pI = p{2}(I,:); pbI = [1:size(pI,1); 1:size(pI,1)]';

%% solve eigenproblem
I = speye(size(D));
tic; [V,L] = eigs(D,M,10,'SM'); toc
[lam,ord] = sort(diag(L),'descend'); V = V(:,ord); lam

%% plot spectrum
figure(1); plot(lam,'*'); axis tight, axis square
xlabel('$k$'); ylabel('$\lambda_k$');

%% plot eigenvector
figure(2), clf; tI = delaunayn(pI); 
w = V(:,3);
plotf(pI,tI,pbI,w/norm(w),1); axis([0 1 0 1]); colorbar
xlabel('$x$'); ylabel('$y$');

%% compute partition
n1 = 200; x1 = linspace(0,1,n1); [X1,Y1] = meshgrid(x1,x1); 
W = eval_p1(p{1},V(:,1:3),[X1(:) Y1(:)]);       % evaluate eigenvectors on grid
idx = kmeans(W, size(W,2));                     % kmeans clustering

%% plot partition
figure(3); clf; surf(X1,Y1,reshape(idx,n1,n1)); view(2); shading flat
axis equal; axis tight; xlabel('$x$'); ylabel('$y$'); colorbar
colormap(70*[ 1 1 1; 2 2 2; 3 3 3]/255);

