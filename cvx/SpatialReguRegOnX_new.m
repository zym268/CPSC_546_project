function beta = SpatialReguRegOnX_new( X, Y, D, lambda, gamma, index)
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

%thre = 0.7;  % threshold on the weight
% first translate the distance matrix to the spatial weight matrix W = 1/D
% W = 1./D;
% W = tril(W,-1);
% Wvalue = sort(W(:),'descend');
% t = min(3000,round(xN*(xN-1)*0.2)); %what does this mean?
% %W(W < Wvalue(t)) = 0;
% W(W<Wvalue(t)) = 0;
% % sparisity representation of the weight
% ind = find(W);
% G2.E = inds2subs(ind, size(W));
% G2.C = W(ind);
% G2.W = G2.C;
% if nargin == 6
%     for i = 1:length(ind)
%         Scale = sqrt(length(index{G2.E(1,i)}).*length(index{G2.E(2,i)}));
%         G2.W(i) = G2.C(i)/Scale;
%     end
% end

% cvx_begin
% toc
% %cvx_quiet(true);
% variable beta(P);
% expression penalty; expression i;  expression j;
% penalty = 0;
% for i = 1:length(ind)
%     penalty = penalty + gamma*G2.W(i)*sum(abs(beta(G2.E(1,i))-beta(G2.E(2,i))));
% end
% minimize( sum(sum(square(Y - X*beta))) + penalty );
% cvx_end
% tic

n = size(X,2);

% difference matrix
e = ones(n,1);
D = spdiags([e -e], 0:1, n,n);

cvx_begin
toc
cvx_solver SDPT3
cvx_quiet false
%cvx_quiet(true);
variable beta(P);
minimize( sum(sum(square(Y - X*beta))) + lambda * sum(abs(D*beta) ));
cvx_end
tic





