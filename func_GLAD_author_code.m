function record = func_GLAD_author_code(n_state, n_pt, n_event, provided_label_mat, loc_cover_mat, event_label)
% test the glad algorithm to run on the loc data

[event_ID, pt_ID] = meshgrid(1:n_event,1:n_pt);
event_ID = event_ID(:);
pt_ID = pt_ID(:);

given_labels = provided_label_mat(:);

% find which labels are actually provided - the loc mat is 1
loc_cover_vec = loc_cover_mat(:);
idx = (loc_cover_vec == 1);

% update
event_ID = event_ID(idx);
pt_ID = pt_ID(idx);
given_labels = given_labels(idx);

% prior
p_z1 =0.5;

% Perform EM to infer pZ, beta, and alpha
[ imageStats, labelerStats ] = em(event_ID, pt_ID, given_labels, p_z1, ones(n_pt,1), ones(n_event,1));

% Infer image labels by thresholding pZ; report accuracy
inferred_labels = imageStats{2} >= 0.5;
disp(sprintf('Accuracy using optimal inference: %f', sum(inferred_labels == event_label.') / n_event));

% % Infer image labels by majority vote; report accuracy
% majorityVoteLabels = zeros(n_event,1);
% for i = 1:n_event
% 	idxs = find(event_ID == i);
% 	majorityVoteLabels(i) = round(sum(given_labels(idxs)) / length(idxs));
% end
% disp(sprintf('Accuracy using majority vote: %f', sum(majorityVoteLabels == event_label.') / n_event));

Z_est = [1-inferred_labels inferred_labels];
Z_est = Z_est.';
% confusion matrix
[~, confmtx] = confmat(n_state, n_event, event_label, [0 1], Z_est.');
[~, ~, ~, pre1, rec1, F, WAF] = binary_f_measure(confmtx);

% estimate a and b
[a_est, b_est] = est_ab_prob(n_pt, Z_est, provided_label_mat, loc_cover_mat);

record.a_est = a_est;
record.b_est = b_est;
record.Z_est = Z_est;
record.confmtx = confmtx;
record.F = F;
record.WAF = WAF;
record.pre = pre1;
record.rec = rec1;