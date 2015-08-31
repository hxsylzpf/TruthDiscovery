function record = func_investment(n_pt, n_event, provided_info_mat, event_label, presence_mat)

n_iter = 100;

g = 1.2;
eps = 1e-4;

% n_max is the max num of different claims for an event

% initialization - belief of claims
result_mv = mv_mul(provided_info_mat, presence_mat);
% unique answers for each claim in each event
uni_ans = result_mv.uni_ans; % a cell (n_event,1)

% trustworthy scores of sources - initialization
T = 0.8*ones(n_pt, 1);

% belief of claims
B = result_mv.ratio_ans; % a cell (n_event,1)

for kk = 1:n_iter
% update trustworthy scores of sources
new_T = zeros(n_pt, 1);

for i = 1:n_pt
    for j = 1:n_event
       % provide a fact
       if presence_mat(i,j) == 1
           % what is the claimed value
           cur_info = provided_info_mat(i,j);
           cur_idx = (cur_info == uni_ans{j});
           cur_belief = B{j}(cur_idx);
           cur_trust_vec = T./sum(presence_mat,2);
           % find who also made this claim
           cur_source = (provided_info_mat(:,j) == cur_info);
           % the sum is only over sources who also made this claim
           new_T(i) = new_T(i) + cur_belief*cur_trust_vec(i)/sum(cur_trust_vec(cur_source));
       end
    end
end

% update belief of claims
new_B = cell(n_event, 1);

for j = 1:n_event
    cur_n_ans = length(B{j});
    for k = 1:cur_n_ans
       % what is the claimed value
       cur_info = uni_ans{j}(k);
       cur_trust_vec = new_T./sum(presence_mat,2);
       % find who also made this claim
       cur_source = (provided_info_mat(:,j) == cur_info);
       new_B{j}(k) = sum(cur_trust_vec(cur_source))^g;
    end
end
    
    if norm(new_T - T) <= eps
        break;
    else
        T = new_T;
        B = new_B;
    end
end

% initiated as 0
est_label = zeros(1, n_event);
est_conf = zeros(1, n_event);
% generate the final label for event
for j = 1:n_event
    cur_belief = new_B{j};
    [est_conf(j), cur_idx] = max(cur_belief);
    est_label(j) = uni_ans{j}(cur_idx);
end

accuracy = sum(est_label == event_label)/length(event_label);
record.accuracy = accuracy;
record.est_label = est_label;
record.est_conf = est_conf;

% % confusion matrix
% [~, confmtx] = confmat(n_state, n_event, event_label, [0 1], est_label)
% [~, ~, ~, pre1, rec1, F, WAF] = binary_f_measure(confmtx)
%    
%     record.confmtx = confmtx;
%     record.F = F;
%     record.WAF = WAF;
%     record.est_label = est_label;
%     record.pre = pre1;
%     record.rec = rec1;
    