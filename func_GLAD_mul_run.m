function [record, max_fval, idx] = func_GLAD_mul_run(n_state, n_pt, n_event, provided_label_mat, event_label, loc_cover_mat, Z_est)

n_run = 5;
max_fval = -inf;
% min_fval = inf;

for i = 1:n_run
    cur_result = func_GLAD(n_state, n_pt, n_event, provided_label_mat, event_label, loc_cover_mat, Z_est);
%     cur_result = func_GLAD_grad_ini(n_state, n_pt, n_event, provided_label_mat, event_label, loc_cover_mat, Z_est);
    cur_result.time = toc;
    
    if cur_result.fval > max_fval % < min_fval
        max_fval = cur_result.fval
        record = cur_result;
        idx = i;
    end
end


