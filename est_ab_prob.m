function [a_est, b_est] = est_ab_prob(n_pt, Z_est, provided_label_mat, loc_cover_mat)

% Z_est: 2*n_event, prob of an event being [0 1]

a_est = zeros(1, n_pt);
b_est = zeros(1, n_pt);

for i = 1:n_pt
    % find which locs have been visited, thus which reports have been provided
    cur_idx = (loc_cover_mat(i,:) == 1);
    cur_Z_est = Z_est(:, cur_idx); % 0 or 1
    % find the reports provided by the ith pt
    cur_provided_label = provided_label_mat(i, cur_idx);
    a_est(i) = cur_provided_label*cur_Z_est(2,:).'/(sum(cur_Z_est(2,:))+1e-10)
    b_est(i) = cur_provided_label*cur_Z_est(1,:).'/(sum(cur_Z_est(1,:))+1e-10)
end