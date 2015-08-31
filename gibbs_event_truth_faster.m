function record = gibbs_event_truth_faster(n_state, n_pt, n_event, provided_label_mat, event_label, Z_est)
% gibbs sampling

% hyperparameter
lam_a1 = 1; lam_a0 = 1;
lam_b1 = 1; lam_b0 = 2;
lam_c1 = 1; lam_c0 = 1e8;
lam_r1 = 2; lam_r0 = 2;
lam_g1 = 2; lam_g0 = 2;

% num of run
n_run = 300;
run_to_use = 240:3:300;

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

%% outer run

for kk = 1:n_run

%% first sample z_j

    %% initial counts for all pts and all events
    pt_count_111 = zeros(n_pt, n_event);
    pt_count_110 = zeros(n_pt, n_event);
    pt_count_101 = zeros(n_pt, n_event);
    pt_count_100 = zeros(n_pt, n_event);
    
        for i = 1:n_pt
            cur_HZX = [H(i,:).' Z X(i,:).'];
            % first count all the numbers
            pt_count_111_all = sum(comp_mat_vec(cur_HZX, [1 1 1])); % a_i
            pt_count_110_all = sum(comp_mat_vec(cur_HZX, [1 1 0])); % 1-a_i
            pt_count_101_all = sum(comp_mat_vec(cur_HZX, [1 0 1])); % b_i
            pt_count_100_all = sum(comp_mat_vec(cur_HZX, [1 0 0])); % 1-b_i
            % then revise to account for the removal of jth event
            cur_idx1 = (Z == 1);
            cur_idx0 = (Z == 0);
            
            pt_count_111(i, :) = pt_count_111_all;
            if pt_count_111_all >= 1
                pt_count_111(i, cur_idx1) = pt_count_111(i, cur_idx1) - 1;
            end
            
            pt_count_110(i, :) = pt_count_110_all;
            if pt_count_110_all >= 1
                pt_count_110(i, cur_idx1) = pt_count_110(i, cur_idx1) - 1;
            end
            
            pt_count_101(i, :) = pt_count_101_all;
            if pt_count_101_all >= 1
                pt_count_101(i, cur_idx0) = pt_count_101(i, cur_idx0) - 1;
            end
            
            pt_count_100(i, :) = pt_count_100_all;
            if pt_count_100_all >= 1
                pt_count_100(i, cur_idx0) = pt_count_100(i, cur_idx0) - 1;
            end
            
        end
        
for j = 1:n_event
    % make a copy of last value
    last_Z = Z;
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
%         cur_HZX = [H(i,:).' Z X(i,:).'];
%         % for counting purpose only; first copy
%         cur_HZX_j = cur_HZX;
%         % remove the jth row - corresponding to the jth event
%         cur_HZX_j(j,:) = [];
% %         cur_HX_j = cur_HZX_j(:, [1 3]);
%         pt_count_111 = sum(comp_mat_vec(cur_HZX_j, [1 1 1])); % a_i
%         pt_count_110 = sum(comp_mat_vec(cur_HZX_j, [1 1 0])); % 1-a_i
%         pt_count_101 = sum(comp_mat_vec(cur_HZX_j, [1 0 1])); % b_i
%         pt_count_100 = sum(comp_mat_vec(cur_HZX_j, [1 0 0])); % 1-b_i
        
        if H(i,j) == 1 && X(i,j) == 1
            cur_p_z0 = (pt_count_101(i,j) + lam_b1)/(pt_count_101(i,j) + pt_count_100(i,j) + lam_b1 + lam_b0);
            cur_p_z1 = (pt_count_111(i,j) + lam_a1)/(pt_count_111(i,j) + pt_count_110(i,j) + lam_a1 + lam_a0);
            
            p_suc_z = p_suc_z*cur_p_z1;
            p_fail_z = p_fail_z*cur_p_z0;
            
        elseif H(i,j) == 1 && X(i,j) == 0
            cur_p_z0 = (pt_count_100(i,j) + lam_b0)/(pt_count_101(i,j) + pt_count_100(i,j) + lam_b1 + lam_b0);
            cur_p_z1 = (pt_count_110(i,j) + lam_a0)/(pt_count_111(i,j) + pt_count_110(i,j) + lam_a1 + lam_a0);
            
            p_suc_z = p_suc_z*cur_p_z1;
            p_fail_z = p_fail_z*cur_p_z0;

        end
            
    end
    
    prop_suc_z = p_suc_z/(p_suc_z + p_fail_z);
    % sample and revise the Z value
    Z(j) = my_bernoulli(1, prop_suc_z); % random('bino', 1, prop_suc_z);
    
    %% once Z(j) is sampled, update all the other counts if Z(j) changes
    cur_idx = setdiff(1:n_event, j);
    if Z(j) == 1 && last_Z(j) == 0 % change from 0 to 1
        pt_count_111(:, cur_idx) = pt_count_111(:,cur_idx) + 1;
        pt_count_110(:, cur_idx) = pt_count_110(:,cur_idx) + 1;
        pt_count_101(:, cur_idx) = max(pt_count_101(:,cur_idx) - 1, zeros(n_pt, n_event-1));
        pt_count_100(:, cur_idx) = max(pt_count_100(:,cur_idx) - 1, zeros(n_pt, n_event-1));
    elseif Z(j) == 0 && last_Z(j) == 1 % change from 1 to 0
        pt_count_111(:, cur_idx) = max(pt_count_111(:,cur_idx) - 1, zeros(n_pt, n_event-1));
        pt_count_110 (:, cur_idx) = max(pt_count_110(:,cur_idx) - 1, zeros(n_pt, n_event-1));
        pt_count_101(:, cur_idx) = pt_count_101(:,cur_idx) + 1;
        pt_count_100(:, cur_idx) = pt_count_100(:,cur_idx) + 1;    
    end
    
end

% record Z
record_Z(:,:,kk) = Z;
    
%% then sample H_ij
for j = 1:n_event
        
    for i = 1:n_pt
        
        % sample H_ij only when X_ij is not provided
        if X(i,j) == 0
        
        %% first part
        % consider the jth column of H
        % first remove the ith element from the jth column
        H_j_i = H(:,j); % jth loc
        H_j_i(i) = []; % how other pt visited this loc
        % count number
        h_count_0 = sum(H_j_i == 0);
        h_count_1 = sum(H_j_i == 1);
        cur_p_h0 = h_count_0 + lam_g0;
        cur_p_h1 = h_count_1 + lam_g1;
   
        % initialization; the value is reset per iteration
        p_suc_h = cur_p_h1;
        p_fail_h = cur_p_h0;
    
        %% second part
        cur_HZX = [H(i,:).' Z X(i,:).'];
        % for counting purpose only; first copy
        cur_HZX_j = cur_HZX;
        % remove the jth row - corresponding to the jth event
        cur_HZX_j(j,:) = [];
        cur_HX_j = cur_HZX_j(:, [1 3]);
        pt_count_111 = sum(comp_mat_vec(cur_HZX_j, [1 1 1])); % a_i
        pt_count_110 = sum(comp_mat_vec(cur_HZX_j, [1 1 0])); % 1-a_i
        pt_count_101 = sum(comp_mat_vec(cur_HZX_j, [1 0 1])); % b_i
        pt_count_100 = sum(comp_mat_vec(cur_HZX_j, [1 0 0])); % 1-b_i
        pt_count_0x1 = sum(comp_mat_vec(cur_HX_j, [0 1])); % c
        pt_count_0x0 = sum(comp_mat_vec(cur_HX_j, [0 0])); % 1-c

        % switch
            if Z(j) == 1 % new_sample
                cur_p_h0 = (pt_count_0x0 + lam_c0)/(pt_count_0x1 + pt_count_0x0 + lam_c1 + lam_c0);
                cur_p_h1 = (pt_count_110 + lam_a0)/(pt_count_111 + pt_count_110 + lam_a1 + lam_a0);

                p_suc_h = p_suc_h*cur_p_h1;
                p_fail_h = p_fail_h*cur_p_h0;

            elseif Z(j) == 0
                cur_p_h0 = (pt_count_0x0 + lam_c0)/(pt_count_0x1 + pt_count_0x0 + lam_c1 + lam_c0);
                cur_p_h1 = (pt_count_100 + lam_b0)/(pt_count_101 + pt_count_100 + lam_b1 + lam_b0);

                p_suc_h = p_suc_h*cur_p_h1;
                p_fail_h = p_fail_h*cur_p_h0;
            end
            
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
    record.record_Z = record_Z;
    record.record_H = record_H;
    record.confmtx = confmtx_gibbs;
    record.F = F_gibbs;
    record.WAF = WAF;
    record.pre = pre1;
    record.rec = rec1;
    record.sample_Z = record_Z(:,:,run_to_use);
    record.sample_H = record_H(:,:,run_to_use);
    
