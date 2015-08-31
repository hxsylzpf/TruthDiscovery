function [h_est, lambda_est, rmse, support] = pt_bias_pre_two_level(X, z, n_l, diff_level_ind)
% calculate pt's bias and precision given z and difficulty level
% diff_level_ind is in [1,2], annotating which z is in level 1 and which in level 2

n_pt = size(X, 1);

h_est = zeros(n_pt, n_l);
lambda_est = zeros(n_pt, n_l);
rmse = zeros(n_pt, n_l);
support = zeros(n_pt, n_l);

for i = 1:n_pt
    for k = 1:n_l
        k
        % provide est for target quantities with diff level as k
        idx_i = (~isnan(X(i,:)) & (diff_level_ind == k));
        if sum(idx_i) > 1
            err = X(i,idx_i) - z(idx_i);
            h_est(i,k) = mean(err);
            cur_var = mean((err - h_est(i,k)).^2);
            lambda_est(i,k) = 1/(cur_var + 1e-2);
            rmse(i,k) = my_rmse(err);
            support(i,k) = sum(idx_i);
        else
            h_est(i,k) = 0;
            lambda_est(i,k) = 1e-3;
            rmse(i,k) = 0;
            support(i,k) = 1e-3;
        end
    end
end

