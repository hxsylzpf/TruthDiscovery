%% model parameter estimation from given data - EM, incremental
%% note: the data instance is provided one by one
% the input is a data table containing [pt_id, fact_id, val]
function record = func_em_ci_mul_pre_bias_ini_incre(data_table, last_val, n_l, ini_z, ini_lambda, ini_h, prior_mu, prior_nu, ....
                  prior_a, prior_b)
% n_l - num of levels
% last_val is a struct, storing all the cached values
% parameter initialization is very important for EM - bad initialization can result in bad estimation

%% max n_pt and n_s
cur_max_n_pt = max(data_table(:,1));
cur_max_n_s = max(data_table(:,2));
    
if ~isempty(last_val)
    n_pt_last = last_val.max_n_pt;
    n_s_last = last_val.max_n_s;
    max_n_pt = max([cur_max_n_pt last_val.max_n_pt]);
    max_n_s = max([cur_max_n_s last_val.max_n_s]);
    n_claim = last_val.n_claim;
    % ini
    n_claim_pt = zeros(max_n_pt,1);
    n_claim_s = zeros(1,max_n_s);
    % copy value
    n_claim_pt(1:n_pt_last) = last_val.n_claim_pt;
    n_claim_s(1:n_s_last) = last_val.n_claim_s;
else
    max_n_pt = cur_max_n_pt;
    max_n_s = cur_max_n_s;
    n_claim_pt = zeros(max_n_pt,1);
    n_claim_s = zeros(1,max_n_s);
end

%% initialization - very important
    s1 = ones(1,n_l)/n_l;
    s2 = zeros(max_n_pt, n_l);
    s3 = zeros(max_n_pt, n_l);
    s4 = zeros(1, max_n_s);
    s5 = zeros(1, max_n_s);
    
    pi_est = ones(1,n_l)/n_l;
    z_est = ini_z;
    lambda_est = ini_lambda; % rand(n_pt, n_l);
    h_est = ini_h; % 0.1*(-0.5+rand(n_pt, n_l));
        
    % copy last val if not empty
    % update happens only for pts and facts in the incremental data
    if ~isempty(last_val)
        
        s1 = last_val.s1;
        s2(1:n_pt_last,:) = last_val.s2;
        s3(1:n_pt_last,:) = last_val.s3;
        s4(1:n_s_last) = last_val.s4;
        s5(1:n_s_last) = last_val.s5;
        
        pi_est = last_val.pi_est;
        h_est(1:n_pt_last,:) = last_val.h_est;
        lambda_est(1:n_pt_last,:) = last_val.lambda_est;
        z_est(1:n_s_last) = last_val.z_est;
    end

%% num of current data insts
n_cur_data = size(data_table, 1);
   
% debug
s1_record = zeros(n_cur_data, n_l);
s2_record = zeros(n_cur_data, max_n_pt);
s3_record = zeros(n_cur_data, max_n_pt);

%% update values based on current data table
% the number of iteration is the same as the number of current data instances
for nn = 1:n_cur_data
    
    %% extract info from data
    cur_pt_id = data_table(nn,1);
    cur_s_id = data_table(nn,2);
    cur_val = data_table(nn,3);
    
    %% update count
    n_claim = n_claim + 1;
    n_claim_pt(cur_pt_id) = n_claim_pt(cur_pt_id) + 1;
    
    % different step sizes
    step = (n_claim)^(-0.6);
    step_pt = (n_claim_pt(cur_pt_id))^(-0.6);
%     step_s = (n_claim_s(cur_s_id))^(-0.6);
    
    %% E-step - calculate expectation and sufficient statistics    
    %% p
    p = zeros(1,n_l);
    for k = 1:n_l
        p(k) = pi_est(k)*my_normpdf(cur_val,z_est(cur_s_id)+h_est(cur_pt_id,k),sqrt(1/lambda_est(cur_pt_id,k)));
    end 
    
    % normalize
    p = p/sum(p);
    
    p = p + 1e-3;
    p = p/sum(p);
        
    %% s1
    s1 = (1-step)*s1 + step*p;
    
    %% s2
    s2(cur_pt_id,:) = (1-step_pt)*s2(cur_pt_id,:) + step_pt*p*cur_val;
    
    %% s3
    s3(cur_pt_id,:) = (1-step_pt)*s3(cur_pt_id,:) + step_pt*p*cur_val^2;
    
    % debug
    s1_record(nn,:) = s1;
    s2_record(nn,:) = s2(:,1).';
    s3_record(nn,:) = s3(:,1).';
    
    %% M-step - estimate model paras - directly reassign
    pi_est = s1;
    
    for iter = 1:100
        h_est(cur_pt_id,:) = (s2(cur_pt_id,:) - s1*z_est(cur_s_id))./s1;
        lambda_est(cur_pt_id,:) = (0.5*s1 + prior_a(cur_pt_id,:) - 1)./ ...
            (0.5*s3(cur_pt_id,:) + 0.5*s1.*(z_est(cur_s_id) + h_est(cur_pt_id,:)) - s2(cur_pt_id,:).*(z_est(cur_s_id) + h_est(cur_pt_id,:))...
            + prior_b(cur_pt_id,:));

        s4(cur_s_id) = s4(cur_s_id) + lambda_est(cur_pt_id,:)*(s2(cur_pt_id,:) - s1.*h_est(cur_pt_id,:)).';
        s5(cur_s_id) = s5(cur_s_id) + lambda_est(cur_pt_id,:)*s1.';

        temp = (prior_mu(cur_s_id)*prior_nu(cur_s_id) + s4(cur_s_id))/(prior_nu(cur_s_id) + s5(cur_s_id));
        
        % check convergence
        if abs(temp - z_est(cur_s_id)) < 1e-3
            break;
        else
            z_est(cur_s_id) = temp;
        end
    end
       
end

%% record the vals after model iteration
    record.pi_est = pi_est;
    record.h_est = h_est;
    record.lambda_est = lambda_est;
    record.z_est = z_est;
       
    record.max_n_pt = max_n_pt;
    record.max_n_s = max_n_s;
    record.n_claim = n_claim;
    record.n_claim_pt = n_claim_pt;
    record.n_claim_s = n_claim_s;
    
    record.s1 = s1;
    record.s2 = s2;
    record.s3 = s3;
    record.s4 = s4;
    record.s5 = s5;
        