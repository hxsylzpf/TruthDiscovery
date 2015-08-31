function record = gibbs_event_truth_mul(n_state, n_pt, n_event, provided_label_mat, event_label, Z_ini, lam_g1, lam_g0)
% gibbs sampling for multinomial truth
% n_state is the number of states
% an event can take a value from 0 to n_state - 1

% hyperparameter
lam_a1 = 5; lam_a0 = 5;
lam_b1 = 5; lam_b0 = 20;
lam_r = 5*ones(n_state,1);

% num of run
n_run = 800;
run_to_use = 600:2:n_run;

% data
X = provided_label_mat;
% initialize H
H = round(rand(n_pt, n_event));
idx = (X > 0);
H(idx) = 1;

% initialize Z
Z = Z_ini; % col vector
% record the samples in order to perform expectation
record_Z = zeros(n_event, 1, n_run);
record_H = zeros(n_pt, n_event, n_run);

%% outer run

for kk = 1:n_run

    kk
%% first sample z_j

for j = 1:n_event
    
    % only used for counting purpose
    Z_j = Z;
    Z_j(j) = [];
    %% first part
    % z_count(i) stores the counts of z = i-1
    % e.g., z_count(1) stores the counts of z = 0
    z_count = zeros(n_state,1);
    for cc = 1:n_state
        z_count(cc) = sum(Z_j == cc-1) + lam_r(cc);
    end
    
    % initialization; the value is reset per iteration
    log_p_z = log(z_count);

    %% second part
    for i = 1:n_pt
        if H(i,j) == 1
            % sampling z = cc - 1
           for cc = 1:n_state
               log_p_z(cc) = log_p_z(cc) + log(phi_mul(i, j, H, Z, X, 1, cc-1, X(i,j), lam_a1, lam_a0, lam_b1, lam_b0, n_state));
           end
        end
    end
    
    prob_z = exp(log_p_z)/sum(exp(log_p_z));
    % sample and revise the Z value
    Z(j) = discreternd(prob_z,1) - 1;
    
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
        cur_p_h0 = h_count_0 + lam_g0(j);
        cur_p_h1 = h_count_1 + lam_g1(j);
   
        % initialization; the value is reset per iteration
        p_suc_h = cur_p_h1;
        p_fail_h = cur_p_h0;
    
        %% second part
            % only update - sample h = 1; as sample h = 0 is only determined by loc pop
            p_suc_h = p_suc_h*phi_mul(i, j, H, Z, X, 1, Z(j), X(i,j), lam_a1, lam_a0, lam_b1, lam_b0, n_state);
            
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
pro_claim = cell(n_event,1);
% only the prob of the considered truth
% e.g., fact1: c1 is correct (in Z_est) p(c1)
% fact2: c3 is correct (in Z_est) p(c3)
pro_truth = zeros(n_event,1);
Z_est = zeros(1, n_event);

    for j = 1:n_event
        cur_Z = record_Z(j,:,:);
        pro_claim{j} = my_count_unique(cur_Z(:), 0:n_state-1);
        pro_claim{j} = pro_claim{j}/sum(pro_claim{j});
        % index
        Z_est(j) = find(pro_claim{j} == max(pro_claim{j}), 1, 'first');
        pro_truth(j) =  pro_claim{j}(Z_est(j));
    end

% turn to label
Z_est = Z_est - 1;
        
% probability of a pt visited a loc
H_est = mean(record_H(:,:,run_to_use), 3);
H_est_round = round(H_est);

%% true positive rate and true negative rate
[tpr, tnr] = mul_truth_tpr_tnr(Z_est, event_label);
acc = sum(Z_est == event_label)/n_event;

%% calculate a and b
a_est_gibbs = zeros(1, n_pt);
b_est_gibbs = zeros(1, n_pt);

for i = 1:n_pt
      
        cur_HZX = [H_est_round(i,:).' Z_est.' X(i,:).'];
    
        pt_count_111 = sum(cur_HZX(:,1)==1 & cur_HZX(:,2)==cur_HZX(:,3) & cur_HZX(:,2)~=0); % a_i
        pt_count_110 = sum(cur_HZX(:,1)==1 & cur_HZX(:,2)~=cur_HZX(:,3) & cur_HZX(:,2)~=0); % 1-a_i
    
        pt_count_101 = sum(cur_HZX(:,1)==1 & cur_HZX(:,2)==0 & cur_HZX(:,3)~=0); % b_i
        pt_count_100 = sum(cur_HZX(:,1)==1 & cur_HZX(:,2)==0 & cur_HZX(:,3)==0); % 1-b_i
    
        a_est_gibbs(i) = (pt_count_111 + lam_a1)/(pt_count_111 + pt_count_110 + lam_a1 + lam_a0);
        b_est_gibbs(i) = (pt_count_101 + lam_b1)/(pt_count_101 + pt_count_100 + lam_b1 + lam_b0);        
end

% err_a_gibbs = mean(abs(a_est_gibbs - a.'));
% err_b_gibbs = mean(abs(b_est_gibbs - b.'));

%% calculate g
g_est_gibbs = zeros(1, n_event);

for j = 1:n_event
    g_est_gibbs(j) = (sum(H_est(:,j)) + lam_g1)/(n_pt + lam_g1 + lam_g0);
end

%% show results

    record.a_est = a_est_gibbs;
    record.b_est = b_est_gibbs;
    record.g_est = g_est_gibbs;
    record.Z_est = Z_est;
    record.H_est = H_est_round;
    record.record_Z = record_Z;
    record.record_H = record_H;
    record.acc = acc;
    record.tpr = tpr;
    record.tnr = tnr;
    record.sample_Z = record_Z(:,:,run_to_use);
    record.sample_H = record_H(:,:,run_to_use);
    
