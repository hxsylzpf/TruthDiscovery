function record = func_gibbs_event_truth_ini_mul_run(n_state, n_pt, n_event, provided_label_mat, event_label, Z_est, n_run)

% hyperparameter
lam_a1 = 1; lam_a0 = 1;
lam_b1 = 1; lam_b0 = 2;
lam_g1 = 2; lam_g0 = 2;

cur_result = gibbs_event_truth(n_state, n_pt, n_event, provided_label_mat, event_label, Z_est);
sample_Z = cur_result.sample_Z;
sample_H = cur_result.sample_H;

for i = 1:n_run
    % new Z_est
    Z_est = rand(1, n_event);
    Z_est = [Z_est; 1-Z_est];
    cur_result = gibbs_event_truth(n_state, n_pt, n_event, provided_label_mat, event_label, Z_est);
    sample_Z = cat(3, sample_Z, cur_result.sample_Z);
    sample_H = cat(3, sample_H, cur_result.sample_H);
end


%% calculate probabilities
% data
X = provided_label_mat;

% probability of an event being true
Z_est1 = mean(sample_Z, 3);
Z_est0 = 1 - Z_est1;
Z_est = [Z_est0 Z_est1];
Z_est1_round = (Z_est1 >= 0.5);
Z_est0_round = 1 - Z_est1_round;
Z_est_round = [Z_est0_round Z_est1_round];

% est_label = (Z_est1 > 0.5);
% acc = 1 - sum(abs(est_label - event_label.'))/n_event

[~, confmtx_gibbs] = confmat(n_state, n_event, event_label, [0 1], Z_est_round);
[~, ~, ~, pre1, rec1, F_gibbs, WAF] = binary_f_measure(confmtx_gibbs);

% probability of a pt visited a loc
H_est1 = mean(sample_H, 3);
H_est1_round = round(H_est1);

%% calculate a and b
a_est_gibbs = zeros(1, n_pt);
b_est_gibbs = zeros(1, n_pt);

for i = 1:n_pt
%         % correspond to pt i only
%         cur_HZX = [H_est1_round(i,:).' Z_est_round(:,2) X(i,:).'];
        
%         % find index
%         La1 = comp_mat_vec(cur_HZX,[1 1 1]);
%         La0 = comp_mat_vec(cur_HZX,[1 1 0]);
%         Lb1 = comp_mat_vec(cur_HZX,[1 0 1]);
%         Lb0 = comp_mat_vec(cur_HZX,[1 0 0]);

% find index
        L1 = (X(i,:).' == 1);
        L0 = (X(i,:).' == 0);

        % weighted sum
        na1 = sum(L1.*Z_est1.*H_est1(i,:).');
        na0 = sum(L0.*Z_est1.*H_est1(i,:).');
        nb1 = sum(L1.*Z_est0.*H_est1(i,:).');
        nb0 = sum(L0.*Z_est0.*H_est1(i,:).');
        
        % turns to probability
        a_est_gibbs(i) = (na1 + lam_a1)/(na1 + na0 + lam_a1+lam_a0);
        b_est_gibbs(i) = (nb1 + lam_b1)/(nb1 + nb0 + lam_b1+lam_b0);        
        
end

% err_a_gibbs = mean(abs(a_est_gibbs - a.'));
% err_b_gibbs = mean(abs(b_est_gibbs - b.'));

%% calculate g
g_est_gibbs = zeros(1, n_event);

for j = 1:n_event
    g_est_gibbs(j) = (sum(H_est1(:,j)) + lam_g1)/(n_pt + lam_g1 + lam_g0);
end

%% show results

idx = find((Z_est1_round - event_label.')~=0)

    record.a_est = a_est_gibbs;
    record.b_est = b_est_gibbs;
    record.g_est = g_est_gibbs;
    record.Z_est = Z_est.';
    record.H_est = H_est1;
    record.confmtx = confmtx_gibbs;
    record.F = F_gibbs;
    record.WAF = WAF;
    record.pre = pre1;
    record.rec = rec1;
    
