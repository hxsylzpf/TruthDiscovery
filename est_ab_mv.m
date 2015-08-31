function [est_a, est_b] = est_ab_mv(n_pt, est_event_label, provided_label_mat, loc_cover_mat)

est_a = zeros(1, n_pt);
est_b = zeros(1, n_pt);

for i = 1:n_pt
    % find which locs have been visited, thus which reports have been provided
    cur_idx = (loc_cover_mat(i,:) == 1);
    cur_event_label = est_event_label(cur_idx); % 0 or 1
    % find the reports provided by the ith pt
    cur_provided_label = provided_label_mat(i, cur_idx);
    idx0 = (cur_event_label == 0);
    idx1 = (cur_event_label == 1);
    est_a(i) = sum(cur_provided_label == 1 & cur_event_label == 1)/(sum(idx1) + 0.01);
    est_b(i) = sum(cur_provided_label == 1 & cur_event_label == 0)/(sum(idx0) + 0.01);
end