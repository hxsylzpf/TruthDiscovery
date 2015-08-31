function [prior_mu, prior_nu] = func_prior_mu_nu(X)

% this function computes prior_mu and prior_nu from data

n_s = size(X,2);

z_est13 = zeros(1,n_s);
% precision
z_pre_est13 = zeros(1,n_s);
% support - num of pts who made a claim on a fact
z_sup_est13 = zeros(1, n_s);
    
for j = 1:n_s
    idx_j = ~isnan(X(:,j));
    z_sup_est13(j) = sum(idx_j);
    % over nonzero eles
    z_est13(j) = median(X(idx_j, j));
    cur_var = mean((X(idx_j,j) - z_est13(j)).^2);
    z_pre_est13(j) = 1/(cur_var + 0.001);
end

prior_mu = z_est13;
prior_nu = z_sup_est13.*z_pre_est13;
