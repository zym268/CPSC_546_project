function [beta,iteration,history] = total_variation(X,Y, lambda, rho, alpha)
% total_variation  Solve total variation minimization via ADMM
%
% [x, history] = total_variation(b, lambda, rho, alpha)
% 
% Solves the following problem via ADMM:
% 
%   minimize  (1/2)||x*beta - Y||_2^2 + lambda * sum_i |beta_{i+1} - beta_i|
%
% where b in R^n.
%
% The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and 
% dual residual norms, and the tolerances for the primal and dual residual 
% norms at each iteration.
% 
% rho is the augmented Lagrangian parameter. 
%
% alpha is the over-relaxation parameter (typical values for alpha are 
% between 1.0 and 1.8).
%
%
% More information can be found in the paper linked at:
% http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%
%This script has been modified for my own problem.

t_start = tic;

%% Global constants and defaults

QUIET    = 1;
MAX_ITER = 2000;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;


%% Data preprocessing

n = size(X,2);

% difference matrix
e = ones(n,1);
D = spdiags([e -e], 0:1, n,n);

%% ADMM solver

beta = zeros(n,1);
z = zeros(n,1);
u = zeros(n,1);

if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

I = speye(n);
DtD = D'*D;

for k = 1:MAX_ITER

    % beta-update
    beta = (X'*X + rho*DtD) \ (X'*Y + rho*D'*(z-u));

    % z-update with relaxation
    zold = z;
    Ax_hat = alpha*D*beta +(1-alpha)*zold;
    z = shrinkage(Ax_hat + u, lambda/rho);

    % u-update
    u = u + Ax_hat - z;

    % diagnostics, reporting, termination checks
    history.objval(k)  = objective(X,Y, lambda, D, beta, z);

    history.r_norm(k)  = norm(D*beta - z);
    history.s_norm(k)  = norm(-rho*D'*(z - zold));
    
    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(D*beta), norm(-z));
    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*D'*u);
    
    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
    end
    
    if (history.r_norm(k) < history.eps_pri(k) && ...
       history.s_norm(k) < history.eps_dual(k))
         break;
    end
end

if ~QUIET
    toc(t_start);
end

iteration=k;

end

function obj = objective(X,Y, lambda, D, beta, z)
    obj = .5*norm(X*beta - Y)^2 + lambda*norm(z,1);
end

function y = shrinkage(a, kappa)
    y = max(0, a-kappa) - max(0, -a-kappa);
end


