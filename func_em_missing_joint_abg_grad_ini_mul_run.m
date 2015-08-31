function [record, min_fval, idx] = func_em_missing_joint_abg_grad_ini_mul_run(n_state, n_pt, n_event, provided_label_mat, event_label, Z_est, a_est, b_est, g_est, n_run)

tic
cur_result = func_em_missing_joint_abg_grad_ini(n_state, n_pt, n_event, provided_label_mat, event_label, Z_est, a_est, b_est, g_est);
cur_result.time = toc;
% keep record of the best one
min_fval = cur_result.fval;
record = cur_result;
idx = 1;

for i = 1:n_run
    % new Z_est
    Z_est = rand(1, n_event);
    Z_est = [Z_est; 1-Z_est];
    tic
    cur_result = func_em_missing_joint_abg_grad_ini(n_state, n_pt, n_event, provided_label_mat, event_label, Z_est, a_est, b_est, g_est);
    cur_result.time = toc;
    
    if cur_result.fval < min_fval
        min_fval = cur_result.fval;
        record = cur_result;
        idx = i+1;
    end
end
