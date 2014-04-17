function beta = SpatialReguRegOnX_SPG( X, Y, D, lambda, gamma)
% using cvx toolbox to perform the spatial regulzaried lasso with G2
% representing the distance weight
% Input:    X ~ N*P 
%           Y ~ N*J
%           D ~ K*K, distance measure, assume to apply the spatial
%           regularization on the first K variables in X
%           lambda  ~sparsity penalty
%           gamma   ~spatial regularied penalty
% Output:   beta learned coefficients  

yN = size(Y,2);  % 
xN = size(D,1); % dimension of spatiall regularied variables
P = size(X,2); % dimension of total variables in X

if nargin < 6,
    for i = 1:xN
        index{i} = i;
    end
end


%SPG algorithm
%initialize, normally set alpha and rho as 1
alpha=1.0;
rho=1.0;
n=size(X,2);
e = ones(n,1);
C = spdiags([e -e], 0:1, n,n);

[beta, obj, density, iter,  time, L]=SPG_linesearch('graph', Y, X, gamma, lambda, C, 'display_iter');