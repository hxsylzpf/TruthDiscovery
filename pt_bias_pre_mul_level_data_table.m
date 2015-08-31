function [h_est, lambda_est, rmse, support] = pt_bias_pre_mul_level_data_table(data_table, z, n_l, diff_level_ind)
% calculate pt's bias and precision given z and difficulty level
% diff_level_ind is in [1,2,3, ...], annotating which z is in level 1, 2, ...

n_pt = max(data_table(:,1));

% note: these are only used for initialization
% the actual returned h and lambda may be much longer
h_est = zeros(n_pt, n_l);
lambda_est = zeros(n_pt, n_l);
rmse = zeros(n_pt, n_l);
support = zeros(n_pt, n_l);

for i = 1:n_pt
    % cur pt
    cur_pt_id = i;
    % row indicator of val provided
    cur_row_idx = (data_table(:,1) == cur_pt_id);
    % s provided by cur pt - integer idx (compressed vec)
    cur_s_by_pt_int_id = data_table(cur_row_idx,2);
    % convert to logical idx with the same length as diff_level_ind
    cur_s_by_pt_idx = zeros(size(diff_level_ind));
    cur_s_by_pt_idx(cur_s_by_pt_int_id) = 1;
    
    % val provided by cur pt - compressed vec
    cur_val = data_table(cur_row_idx,3);
    
    % val provided - a vector of the same size as z
    cur_val_vec = zeros(size(z));
    cur_val_vec(cur_s_by_pt_int_id) = cur_val.';
        
    for k = 1:n_l
        k
        % provide est for target quantities with diff level k
        idx_i = ((cur_s_by_pt_idx) & (diff_level_ind == k));
        if sum(idx_i) > 3
            err = cur_val_vec(idx_i) - z(idx_i);
            h_est(cur_pt_id,k) = mean(err);
            cur_var = mean((err - h_est(cur_pt_id,k)).^2);
            lambda_est(cur_pt_id,k) = 1/(cur_var + 1e-2);
            rmse(cur_pt_id,k) = my_rmse(err);
            support(cur_pt_id,k) = sum(idx_i);
        else
            h_est(cur_pt_id,k) = 0;
            lambda_est(cur_pt_id,k) = 1e-6;
            rmse(cur_pt_id,k) = 0;
            support(cur_pt_id,k) = 1e-6;
        end
    end
end

