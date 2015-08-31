function record = gibbs_event_truth_bi_online(n_state, n_pt, n_event, provided_label_mat, event_label, s_est, g_est, a_est, b_est)

% TSE online

%% post prob
prob = zeros(n_state, n_event);

for j = 1:n_event
    prob(1,j) = 1-s_est; % 0
    prob(2,j) = s_est; % 1
    for i = 1:n_pt
        cur_claim = provided_label_mat(i,j);
        % p(z=0)
        prob(1,j) = prob(1,j)*(g_est(j)*b_est(i)^(cur_claim)*(1-b_est(i))^(1-cur_claim) ...
                             +(1-g_est(j))*0^(cur_claim)*1^(1-cur_claim));
        % p(z=1)
        prob(2,j) = prob(2,j)*(g_est(j)*a_est(i)^(cur_claim)*(1-a_est(i))^(1-cur_claim) ...
                             +(1-g_est(j))*0^(cur_claim)*1^(1-cur_claim));               
    end
end

prob = prob_mat_nlz(prob, 'col');

Z_est = round(prob(2,:));

[~, confmtx_gibbs] = confmat(n_state, n_event, event_label, [0 1], Z_est);
[~, ~, ~, pre1, rec1, F_gibbs, WAF] = binary_f_measure(confmtx_gibbs);

record.Z_est = Z_est.';
record.confmtx = confmtx_gibbs;
record.F = F_gibbs;
record.WAF = WAF;
record.pre = pre1;
record.rec = rec1;
