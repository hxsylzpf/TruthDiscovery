function record = gibbs_truth_fast(n_state, n_pt, n_event, provided_label_mat, loc_cover_mat, event_label, Z_est)
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

% global counts of Z = 0 and Z = 1
n_z0 = sum(1 - Z);
n_z1 = sum(Z);

% global counts of zx
n_zx_00 = zeros(n_pt, 1);
n_zx_01 = zeros(n_pt, 1);
n_zx_10 = zeros(n_pt, 1);
n_zx_11 = zeros(n_pt, 1);

for i = 1:n_pt
    n_zx_00(i) = sum((Z.' == 0) & (X(i,:) == 0));
    n_zx_01(i) = sum((Z.' == 0) & (X(i,:) == 1));
    n_zx_10(i) = sum((Z.' == 1) & (X(i,:) == 0));
    n_zx_11(i) = sum((Z.' == 1) & (X(i,:) == 1));
end

% record the samples in order to perform expectation
record_Z = zeros(n_event, 1, n_run);

%% outer run

for kk = 1:n_run

    kk
    
%% first sample z_j

for j = 1:n_event
    
    %% first part

    % before sampling
    last_z = Z(j);
    
    % step 1: calculate local counts
    if Z(j) == 0
       n_z0_local = n_z0 - 1;
       n_z1_local = n_z1;
    elseif Z(j) == 1
       n_z0_local = n_z0;
       n_z1_local = n_z1 - 1;
    end
    
    cur_p_z0 = n_z0_local + lam_r0;
    cur_p_z1 = n_z1_local + lam_r1;
    
    % initialization; the value is reset per iteration
    p_suc_z = cur_p_z1;
    p_fail_z = cur_p_z0;
    
    %% second part
    for i = 1:n_pt
        
        % update only when a report is provided
        if loc_cover_mat(i,j) == 1

           if Z(j) == 0 && X(i,j) == 0
               n_zx_00_local = n_zx_00(i) - 1;
               n_zx_01_local = n_zx_01(i);
               n_zx_10_local = n_zx_10(i);
               n_zx_11_local = n_zx_11(i);
           elseif Z(j) == 0 && X(i,j) == 1
               n_zx_00_local = n_zx_00(i);
               n_zx_01_local = n_zx_01(i) - 1;
               n_zx_10_local = n_zx_10(i);
               n_zx_11_local = n_zx_11(i);
           elseif Z(j) == 1 && X(i,j) == 0
               n_zx_00_local = n_zx_00(i);
               n_zx_01_local = n_zx_01(i);
               n_zx_10_local = n_zx_10(i) - 1;
               n_zx_11_local = n_zx_11(i);
           elseif Z(j) == 1 && X(i,j) == 1
               n_zx_00_local = n_zx_00(i);
               n_zx_01_local = n_zx_01(i);
               n_zx_10_local = n_zx_10(i);
               n_zx_11_local = n_zx_11(i) - 1;           
           end
             
            if X(i,j) == 1
                cur_p_z0 = (n_zx_01_local + lam_b1)/(n_zx_01_local + n_zx_00_local + lam_b1 + lam_b0);
                cur_p_z1 = (n_zx_11_local + lam_a1)/(n_zx_11_local + n_zx_10_local + lam_a1 + lam_a0);

                p_suc_z = p_suc_z*cur_p_z1;
                p_fail_z = p_fail_z*cur_p_z0;

            elseif X(i,j) == 0
                cur_p_z0 = (n_zx_00_local + lam_b0)/(n_zx_01_local + n_zx_00_local + lam_b1 + lam_b0);
                cur_p_z1 = (n_zx_10_local + lam_a0)/(n_zx_11_local + n_zx_10_local + lam_a1 + lam_a0);

                p_suc_z = p_suc_z*cur_p_z1;
                p_fail_z = p_fail_z*cur_p_z0;

            end
        end
    end
    
    prop_suc_z = p_suc_z/(p_suc_z + p_fail_z);
    % sample and revise the Z value
    cur_z = my_bernoulli(1, prop_suc_z); % random('bino', 1, prop_suc_z);
    Z(j) = cur_z;
    
    %% update all the counts
    [n_z0, n_z1] = func_update_z_count(n_z0, n_z1, last_z, cur_z);
    
    for i = 1:n_pt
        [n_zx_00, n_zx_01, n_zx_10, n_zx_11] = func_update_zx_count(i, n_zx_00, n_zx_01, n_zx_10, n_zx_11, last_z, cur_z, X(i,j));
    end
    
    
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
