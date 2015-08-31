function [h_est, lambda_est, rmse, support] = pt_bias_pre_data_obj(data_obj, z)
% calculate pt's bias and precision given z

n_pt = data_obj.n_pt;

h_est = zeros(n_pt, 1);
lambda_est = zeros(n_pt, 1);
rmse = zeros(n_pt, 1);
% total num of claims a pt made
support = zeros(n_pt,1);

for i = 1:n_pt
    cur_x = data_obj.X_res(i,:);
    idx_i = (~isnan(cur_x));
    err = cur_x(idx_i) - z(idx_i);
    h_est(i) = mean(err);
    cur_var = mean((err - h_est(i)).^2);
    lambda_est(i) = 1/(cur_var + 1e-2);
    rmse(i) = my_rmse(err);
    support(i) = sum(idx_i);
end
