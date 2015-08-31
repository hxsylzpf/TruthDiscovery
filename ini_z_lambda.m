function [ini_z_mu, ini_z_nu, ini_lambda] = ini_z_lambda(X)

[n_pt, n_s] = size(X);

ini_z_mu = zeros(1,n_s);
ini_z_nu = zeros(1,n_s);
ini_lambda = zeros(n_pt, 1);

for j = 1:n_s
    idx_j = (~isnan(X(:,j)));
    % over nonzero eles
    ini_z_mu(j) = median(X(idx_j, j));
    ini_z_nu(j) = 1/(var(X(idx_j, j)) + 0.001);
end

% incorrect
for i = 1:n_pt
    idx_i = (~isnan(X(i,:)));
    % over nonzero eles
    ini_lambda(i) = 1/(var(X(i,idx_i) - ini_z_mu(idx_i)) + 0.001);
end