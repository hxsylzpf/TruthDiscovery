function [z_est, pre_est] = tf_median(X)

n_s = size(X,2);
z_est = zeros(1,n_s);
pre_est = zeros(1,n_s);

for j = 1:n_s
    idx_j = ~isnan(X(:,j));
    z_est(j) = median(X(idx_j, j));
    cur_var = mean((X(idx_j,j) - z_est(j)).^2);
    pre_est(j) = 1/(cur_var + 1e-4);
end