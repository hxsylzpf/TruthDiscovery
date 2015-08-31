function record = gibbs_event_truth_psn_online(n_state, n_pt, n_event, provided_label_mat, event_label, s_est, a_est, b_est, temp_loc_cover_mat)

% PTSE online

%% infer f
n_pt_temp = size(temp_loc_cover_mat,1);

sim_mat = zeros(n_pt, n_pt_temp);

for i = 1:n_pt
    for k = 1:n_pt_temp
        temp1 = provided_label_mat(i,:);
        temp2 = temp_loc_cover_mat(k,:);
        val1 = sum(temp1 & temp2);
        val2 = sum(temp1 | temp2);
        sim_mat(i,k) = val1/val2;
    end
end

tic
% predict
psn_loc_visit_tend = zeros(n_pt, n_event);

for i = 1:n_pt
    for j = 1:n_event
        psn_loc_visit_tend(i,j) = sim_mat(i,:)*temp_loc_cover_mat(:,j)/sum(sim_mat(i,:));
    end
end

%% post prob
prob = zeros(n_state, n_event);

for j = 1:n_event
    prob(1,j) = 1-s_est; % 0
    prob(2,j) = s_est; % 1
    for i = 1:n_pt
        cur_claim = provided_label_mat(i,j);
        % p(z=0)
        prob(1,j) = prob(1,j)*(psn_loc_visit_tend(i,j)*b_est(i)^(cur_claim)*(1-b_est(i))^(1-cur_claim) ...
                             +(1-psn_loc_visit_tend(i,j))*0^(cur_claim)*1^(1-cur_claim));
        % p(z=1)
        prob(2,j) = prob(2,j)*(psn_loc_visit_tend(i,j)*a_est(i)^(cur_claim)*(1-a_est(i))^(1-cur_claim) ...
                             +(1-psn_loc_visit_tend(i,j))*0^(cur_claim)*1^(1-cur_claim));               
    end
end

prob = prob_mat_nlz(prob, 'col');

Z_est = round(prob(2,:));
cpu_time = toc;

[~, confmtx_gibbs] = confmat(n_state, n_event, event_label, [0 1], Z_est);
[~, ~, ~, pre1, rec1, F_gibbs, WAF] = binary_f_measure(confmtx_gibbs);

record.Z_est = Z_est.';
record.confmtx = confmtx_gibbs;
record.F = F_gibbs;
record.WAF = WAF;
record.pre = pre1;
record.rec = rec1;
record.cpu_time = cpu_time;

