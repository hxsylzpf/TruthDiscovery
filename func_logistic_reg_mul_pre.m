function [est, norm_temp_est] = func_logistic_reg_mul_pre(n_l, w, w0, new_data)

n_data = size(new_data,1);

temp_est = zeros(n_data, n_l);

for k = 1:n_l
    temp_est(:,k) = exp(w(k)*new_data+w0(k));
end

norm_temp_est = prob_mat_nlz(temp_est, 'row');

[~, est] = max(norm_temp_est, [], 2);
