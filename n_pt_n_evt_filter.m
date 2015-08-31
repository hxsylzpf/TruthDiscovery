function [provided_label_mat_f, loc_cover_mat_f, report_loc_f, idx_evt_w_enough_pt] = n_pt_n_evt_filter(min_n_pt, min_n_evt, provided_label_mat, loc_cover_mat, report_loc)

% an event should have min_n_pt report
n_pt_per_evt = sum(provided_label_mat,1);

idx_evt_w_enough_pt = (n_pt_per_evt >= min_n_pt);

report_loc_f = report_loc(idx_evt_w_enough_pt,:);

%% update provided_label_mat - col evts
provided_label_mat_f = provided_label_mat(:, idx_evt_w_enough_pt);
loc_cover_mat_f = loc_cover_mat(:, idx_evt_w_enough_pt);

%% filter out pts with too few reports
idx_pt_w_enough_report = sum(provided_label_mat_f,2) >= min_n_evt;
%% update provided_label_mat - row pts
provided_label_mat_f = provided_label_mat_f(idx_pt_w_enough_report, :);
loc_cover_mat_f = loc_cover_mat_f(idx_pt_w_enough_report, :);

