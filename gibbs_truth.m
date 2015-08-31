function record = gibbs_truth(n_state, n_pt, n_event, provided_label_mat, loc_cover_mat, event_label, Z_est)
% gibbs sampling

% hyperparameter
lam_a1 = 5; lam_a0 = 5;
lam_b1 = 5; lam_b0 = 20;
lam_r1 = 5; lam_r0 = 5;

% num of run
n_run = 400;
run_to_use = 300:2:n_run;

% data
X = provided_label_mat;

% initialize Z
% Z = round(rand(n_event, 1));
Z = round(Z_est(2,:).');
% record the samples in order to perform expectation
record_Z = zeros(n_event, 1, n_run);

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
    p_suc_z = cur_p_z1;
    p_fail_z = cur_p_z0;
    
    %% second part
    for i = 1:n_pt
        
        % update only when report is provided
        if loc_cover_mat(i,j) == 1
            cur_ZX = [Z X(i,:).'];
            % remove the jth row - corresponding to the jth event
            cur_ZX(j,:) = [];
    %         cur_HX_j = cur_HZX_j(:, [1 3]);
%             pt_count_11 = sum(comp_mat_vec(cur_ZX, [1 1])); % a_i
%             pt_count_10 = sum(comp_mat_vec(cur_ZX, [1 0])); % 1-a_i
%             pt_count_01 = sum(comp_mat_vec(cur_ZX, [0 1])); % b_i
%             pt_count_00 = sum(comp_mat_vec(cur_ZX, [0 0])); % 1-b_i
 
            pt_count_11 = sum(cur_ZX(:,1)==1 & cur_ZX(:,2)==1); % a_i
            pt_count_10 = sum(cur_ZX(:,1)==1 & cur_ZX(:,2)==0); % 1-a_i
            pt_count_01 = sum(cur_ZX(:,1)==0 & cur_ZX(:,2)==1); % b_i
            pt_count_00 = sum(cur_ZX(:,1)==0 & cur_ZX(:,2)==0); % 1-b_i
            
            if X(i,j) == 1
                cur_p_z0 = (pt_count_01 + lam_b1)/(pt_count_01 + pt_count_00 + lam_b1 + lam_b0);
                cur_p_z1 = (pt_count_11 + lam_a1)/(pt_count_11 + pt_count_10 + lam_a1 + lam_a0);

                p_suc_z = p_suc_z*cur_p_z1;
                p_fail_z = p_fail_z*cur_p_z0;

            elseif X(i,j) == 0
                cur_p_z0 = (pt_count_00 + lam_b0)/(pt_count_01 + pt_count_00 + lam_b1 + lam_b0);
                cur_p_z1 = (pt_count_10 + lam_a0)/(pt_count_11 + pt_count_10 + lam_a1 + lam_a0);

                p_suc_z = p_suc_z*cur_p_z1;
                p_fail_z = p_fail_z*cur_p_z0;

            end
        end
    end
    
    prop_suc_z = p_suc_z/(p_suc_z + p_fail_z);
    % sample and revise the Z value
    Z(j) = my_bernoulli(1, prop_suc_z); % random('bino', 1, prop_suc_z);
    
end

% record Z
record_Z(:,:,kk) = Z;

% outer iter
end

%% calculate probabilities
% probability of an event being true
Z_est1 = mean(record_Z(:,:,run_to_use), 3);
Z_est0 = 1 - Z_est1;
Z_est = [Z_est0 Z_est1];
Z_est1_round = (Z_est1 >= 0.5);
Z_est0_round = 1 - Z_est1_round;

% est_label = (Z_est1 > 0.5);
% acc = 1 - sum(abs(est_label - event_label.'))/n_event
Z_est_round = [Z_est0_round Z_est1_round];

[~, confmtx_gibbs] = confmat(n_state, n_event, event_label, [0 1], Z_est_round);
[~, ~, ~, pre1, rec1, F_gibbs, WAF] = binary_f_measure(confmtx_gibbs);


%% calculate a and b
a_est_gibbs = zeros(1, n_pt);
b_est_gibbs = zeros(1, n_pt);

for i = 1:n_pt
%         % correspond to pt i only
%         cur_ZX = [Z_est_round(:,2) X(i,:).'];
%         
%         % find index
%         La1 = comp_mat_vec(cur_ZX,[1 1]);
%         La0 = comp_mat_vec(cur_ZX,[1 0]);
%         Lb1 = comp_mat_vec(cur_ZX,[0 1]);
%         Lb0 = comp_mat_vec(cur_ZX,[0 0]);
        
% find index
        L1 = (X(i,:).' == 1);
        L0 = (X(i,:).' == 0);

        % weighted sum
        na1 = sum(L1.*Z_est1);
        na0 = sum(L0.*Z_est1);
        nb1 = sum(L1.*Z_est0);
        nb0 = sum(L0.*Z_est0);
        
        % turns to probability
        a_est_gibbs(i) = (na1 + lam_a1)/(na1 + na0 + lam_a1+lam_a0);
        b_est_gibbs(i) = (nb1 + lam_b1)/(nb1 + nb0 + lam_b1+lam_b0);        
        
end

%% show results

idx = find((round(Z_est1) - event_label.')~=0)

    record.a_est = a_est_gibbs;
    record.b_est = b_est_gibbs;
    record.Z_est = Z_est;
    record.confmtx = confmtx_gibbs;
    record.F = F_gibbs;
    record.WAF = WAF;
    record.pre = pre1;
    record.rec = rec1;
