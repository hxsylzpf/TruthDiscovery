function record = gibbs_event_truth_psn(n_state, n_pt, n_event, provided_label_mat, event_label, Z_est, temp_loc_cover_mat)
% gibbs sampling for binary truth
% it is the same as gibbs_event_truth
% except it is rewritten in terms of a phi function
% temp_loc_cover_mat is a template loc_cover_mat from a set of template pts

% hyperparameter
lam_a1 = 5; lam_a0 = 5;
lam_b1 = 5; lam_b0 = 20;
lam_r1 = 5; lam_r0 = 5;

% num of run
n_run = 400;
run_to_use = 300:2:n_run;

% data
X = provided_label_mat;
% initialize H
H = round(rand(n_pt, n_event));
idx = (X == 1);
H(idx) = 1;

% initialize Z
% Z = round(rand(n_event, 1));
Z = round(Z_est(2,:).');
% Z = ones(n_event, 1);
% record the samples in order to perform expectation
record_Z = zeros(n_event, 1, n_run);
record_H = zeros(n_pt, n_event, n_run);

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

% predict
psn_loc_visit_tend = zeros(n_pt, n_event);

for i = 1:n_pt
    for j = 1:n_event
        if provided_label_mat(i,j) == 0;
           psn_loc_visit_tend(i,j) = sim_mat(i,:)*temp_loc_cover_mat(:,j)/sum(sim_mat(i,:));
        end
    end
end

%% outer run

for kk = 1:n_run

    kk
%% first sample z_j

for j = 1:n_event
    
    % only used for counting purpose
    Z_j = Z;
    Z_j(j) = [];
    %% first part
    z_count_0 = sum(Z_j == 0);
    z_count_1 = sum(Z_j == 1);
    cur_p_z0 = z_count_0 + lam_r0;
    cur_p_z1 = z_count_1 + lam_r1;
    
    % initialization; the value is reset per iteration
    log_p_suc_z = log(cur_p_z1);
    log_p_fail_z = log(cur_p_z0);

    %% second part
    for i = 1:n_pt
        if H(i,j) == 1
           % sample z = 1
           log_p_suc_z = log_p_suc_z + log(phi_bi(i, j, H, Z, X, 1, 1, X(i,j),lam_a1, lam_a0, lam_b1, lam_b0));
           % sample z = 0
           log_p_fail_z = log_p_fail_z + log(phi_bi(i, j, H, Z, X, 1, 0, X(i,j), lam_a1, lam_a0, lam_b1, lam_b0));
        end
    end
    
    p_suc_z = exp(log_p_suc_z);
    p_fail_z = exp(log_p_fail_z);
    prop_suc_z = p_suc_z/(p_suc_z + p_fail_z);
    % sample and revise the Z value
    Z(j) = my_bernoulli(1, prop_suc_z); % random('bino', 1, prop_suc_z);
    
end

% record Z
record_Z(:,:,kk) = Z;
    
%% then sample H_ij
for j = 1:n_event
        
    for i = 1:n_pt
        
        % sample H_ij only when X_ij is not provided
        if X(i,j) == 0
        
            % only update - sample h = 1; as sample h = 0 is only determined by loc pop
            p_suc_h = psn_loc_visit_tend(i,j) * phi_bi(i, j, H, Z, X, 1, Z(j), 0, lam_a1, lam_a0, lam_b1, lam_b0);
            p_fail_h = 1 - psn_loc_visit_tend(i,j);
            
            % sample a new H_ij
            prop_suc_h = p_suc_h/(p_suc_h + p_fail_h);
                        
            H(i,j) = my_bernoulli(1, prop_suc_h); % random('bino', 1, prop_suc_h);    
        end
        
    end
    
end

% record H
record_H(:,:,kk) = H;

% if sum(abs(Z - last_Z)) + sum(sum(abs(last_H - H))) == 0
%     break;
% end
% %% record a copy of current values
% last_Z = Z;
% last_H = H;

% outer iter
end

%% calculate probabilities
% probability of an event being true
Z_est1 = mean(record_Z(:,:,run_to_use), 3);
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
H_est1 = mean(record_H(:,:,run_to_use), 3);
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

%% show results

idx = find((Z_est1_round - event_label.')~=0)

    record.a_est = a_est_gibbs;
    record.b_est = b_est_gibbs;
    record.s_est = sum(Z_est1_round == 1)/n_event;
    record.Z_est = Z_est.';
    record.H_est = H_est1;
    record.record_Z = record_Z;
    record.record_H = record_H;
    record.confmtx = confmtx_gibbs;
    record.F = F_gibbs;
    record.WAF = WAF;
    record.pre = pre1;
    record.rec = rec1;
    record.sample_Z = record_Z(:,:,run_to_use);
    record.sample_H = record_H(:,:,run_to_use);
    
