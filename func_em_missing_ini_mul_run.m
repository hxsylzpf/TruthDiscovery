function [record, max_fval, idx, fval_record, F_record, all_record] = func_em_missing_ini_mul_run(n_state, n_pt, n_event, provided_label_mat, event_label, loc_cover_mat, Z_est)

n_run = 0;

fval_record = zeros(n_run + 1,1);
F_record = zeros(n_run + 1,1);
all_record = struct;

tic
cur_result = func_em_missing_ini(n_state, n_pt, n_event, provided_label_mat, event_label, loc_cover_mat, Z_est);
cur_result.time = toc;
% keep record of the best one
max_fval = cur_result.fval;
record = cur_result;
idx = 1;

fval_record(1) = cur_result.fval;
F_record(1) = cur_result.F;
all_record(1).result = cur_result;

for i = 1:n_run
    % new Z_est
    rng(i);
    Z_est = rand(1, n_event);
    Z_est = [Z_est; 1-Z_est];
    tic
    cur_result = func_em_missing_ini(n_state, n_pt, n_event, provided_label_mat, event_label, loc_cover_mat, Z_est);
    cur_result.time = toc;
    
    fval_record(i+1) = cur_result.fval;
    F_record(i+1) = cur_result.F;
    
    if cur_result.fval > max_fval
        max_fval = cur_result.fval;
        record = cur_result;
        idx = i+1;
    end
end

% % if negative happened; use the first - with provided ini
% if record.F < 0.5
%     record = all_record(1).result;
% end