function [lambda_est, rmse] = pt_pre(X, z)
% calculate pt's bias and precision given z

n_pt = size(X, 1);

lambda_est = zeros(n_pt, 1);
rmse = zeros(n_pt, 1);

for i = 1:n_pt
    idx_i = (~isnan(X(i,:)));
    err = X(i,idx_i) - z(idx_i);
    cur_var = mean(err.^2);
    lambda_est(i) = 1/(cur_var + 1e-2);
    rmse(i) = my_rmse(err);
end