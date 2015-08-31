function [record, min_fval, idx, fval_record, F_record] = func_em_missing_joint_abw_grad_event_ini_mul_run(n_state, n_pt, n_event, n_feature, provided_label_mat,  feature_mat, event_label, Z_est, a_est, b_est, w_est)

% n_run is not needed; automatically provided

w_c = 1./(7:-2:3);

% lambda = n_pt;

n_run = length(w_c);

fval_record = zeros(n_run + 1,1);
F_record = zeros(n_run + 1,1);

% the first run uses the para provided in the main file
tic
cur_result = func_em_missing_joint_abw_grad_event_ini(n_state, n_pt, n_event, n_feature, provided_label_mat,  feature_mat, event_label, Z_est, a_est, b_est, w_est);
cur_result.time = toc;
% cur_w = cur_result.w_est
% norm(cur_w)
% keep record of the best one
min_fval = cur_result.fval;
record = cur_result;
idx = 1;

    fval_record(1) = cur_result.fval;
    F_record(1) = cur_result.F;
%     cur_w = cur_result.w_est;
%     cur_w = reshape(cur_w, n_feature, n_pt);
    
%     for ii = 1:n_pt
% %    cur_w = x(2*n_pt + (ii-1)*n_feature + 1 : 2*n_pt + ii*n_feature);
%     fval_record(1) = fval_record(1) - lambda*norm(cur_w(:,ii), 1)
%     end  
    
for i = 1:n_run
    % new w_est
    W_est = w_c(i)*ones(n_feature, n_event);
    w_est = W_est(:).';
    tic
    cur_result = func_em_missing_joint_abw_grad_event_ini(n_state, n_pt, n_event, n_feature, provided_label_mat,  feature_mat, event_label, Z_est, a_est, b_est, w_est);
    cur_result.time = toc;
%     cur_w = cur_result.w_est
%     norm(cur_w)
%     cur_w = reshape(cur_w, n_feature, n_pt)
    fval_record(i+1) = cur_result.fval;
    F_record(i+1) = cur_result.F;
    
%     
%     for ii = 1:n_pt
% %    cur_w = x(2*n_pt + (ii-1)*n_feature + 1 : 2*n_pt + ii*n_feature);
%     fval_record(i+1) = fval_record(i+1) - lambda*norm(cur_w(:,ii), 1)
%     end  
    
%     if cur_result.fval < min_fval
    if fval_record(i+1) < min_fval
        min_fval = cur_result.fval;
        record = cur_result;
        idx = i+1;
    end
end


