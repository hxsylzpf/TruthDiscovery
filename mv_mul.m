function result = mv_mul(provided_info_mat, presence_mat)
% this func is majority voting on arbitrary possible answers on arbitraty events
% presence_mat must be 1 or 0

[n_pt, n_evt] = size(provided_info_mat);

uni_ans = cell(n_evt, 1);
vote_ans = cell(n_evt, 1);
ratio_ans = cell(n_evt, 1);
predict_val = zeros(n_evt, 1);
predict_conf = zeros(n_evt, 1);
% record whether a pt's answer is the same as majority voting
precision_mat = zeros(n_pt, n_evt);
precision_rate = zeros(n_pt, 1);

for j = 1:n_evt
   % find unique among all the provided values
   idx_j = (presence_mat(:,j) == 1);
   uni_ans{j} = unique(provided_info_mat(idx_j,j));
   vote_ans{j} = mycount_unique(provided_info_mat(idx_j,j), uni_ans{j});
   ratio_ans{j} = vote_ans{j}/sum(vote_ans{j});
   [~, idx] = max(ratio_ans{j});
   predict_val(j) = uni_ans{j}(idx);
   predict_conf(j) = ratio_ans{j}(idx);
   precision_mat(:,j) = ((provided_info_mat(:,j) == predict_val(j))) & (presence_mat(:,j) == 1);
end

for i = 1:n_pt
    precision_rate(i) = sum(precision_mat(i,:))/sum(presence_mat(i,:));
end

result.predict_val = predict_val;
result.predict_conf = predict_conf;
result.uni_ans = uni_ans;
result.ratio_ans = ratio_ans;
result.precision_mat = precision_mat;
result.precision_rate = precision_rate;
