%% model parameter estimation from given data - EM, iterative
% the input is a data table containing [pt_id, fact_id, val]
function record = func_em_ci_mul_pre_bias_ini_iter(data_table, last_val, n_l, ini_z, ini_lambda, ini_h, prior_mu, prior_nu, ....
                  prior_a, prior_b)
% n_pt_sf - total num of pts so far - this value increases over time
% n_s_sf - total num of facts so far - this value increases over time
% n_l - num of levels
% last_val is a struct, storing all the cached values
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% initialization
n_iter = 100;
pi_est = 1/n_l*ones(1, n_l);
outer_n_iter = 0;

% log_val = -inf;
% log_val_record =[];
inner_n_iter_record = [];

% initialize para values
z_est = ini_z;
% [z_est, ~, ~] = ini_z_lambda(X);
lambda_est = ini_lambda; % rand(n_pt, n_l);
h_est = ini_h; % 0.1*(-0.5+rand(n_pt, n_l));

%% max n_pt and n_s
cur_uni_pt = unique(data_table(:,1));
cur_uni_s = unique(data_table(:,2));

if ~isempty(last_val)
%     n_pt_last = last_val.max_n_pt;
    n_s_last = last_val.max_n_s;
    max_n_pt = max([cur_uni_pt; last_val.max_n_pt]);
    max_n_s = max([cur_uni_s; last_val.max_n_s]);
else
    max_n_pt = max(cur_uni_pt);
    max_n_s = max(cur_uni_s);
end

%% num of current data insts
n_cur_data = size(data_table, 1);
   
%% update values based on current data table    
for kk = 1:n_iter
    
    %% E-step - calculate expectation

    %% alpha
    % use log to calculate to avoid num underflow
    
    % allocate space
    log_alpha = zeros(max_n_s, n_l);
    
    % copy last val
    if ~isempty(last_val)
        log_alpha(1:n_s_last,:) = last_val.log_alpha;
    end
    
    % read data table rows one by one
    for nn = 1:n_cur_data
        cur_pt_id = data_table(nn,1);
        cur_s_id = data_table(nn,2);
        cur_val = data_table(nn,3);
        for k = 1:n_l              
            log_alpha(cur_s_id,k) = log_alpha(cur_s_id,k) + log(my_normpdf(cur_val,z_est(cur_s_id)+h_est(cur_pt_id,k),sqrt(1/lambda_est(cur_pt_id,k))));
        end
    end
        
    % only process current n_uni_cur_s
    % convert back to num, also prevent numerical error
    alpha = exp(log_alpha);       
    
    % normalize
    gamma = prob_mat_nlz(repmat(pi_est,max_n_s,1).*alpha, 'row');
    
%     pi_est
%     alpha
%     sum(alpha)

    %% M-step - estimate model paras
    new_pi_est = (sum(gamma)+1)/(max_n_s+n_l);
    
    %% the last updated z value is used as the initial values for para reestimation
    inner_record = func_solve_z_h_lambda_iter(data_table, last_val, n_l, gamma, z_est, prior_mu, prior_nu, prior_a, prior_b);
    new_z_est = inner_record.z_est;
    new_h_est = inner_record.h_est;
    new_lambda_est = inner_record.lambda_est;
%     inner_n_iter_record = [inner_n_iter_record record.n_iter_to_converge];
   
%     new_log_val = func_ci_multiple_Q_bias(n_pt, n_s, n_l, X, new_pi_est, gamma, new_h_est, new_lambda_est, new_z_est, ...
%         prior_mu, prior_nu, prior_mu_h, prior_nu_h, prior_a, prior_b);

%     new_log_val = func_ci_multiple_bias_incmpl_logll(n_pt, n_s, n_l, X, new_pi_est, new_h_est, new_lambda_est, new_z_est, ...
%         prior_mu, prior_nu, prior_mu_h, prior_nu_h, prior_a, prior_b);

%     log_val_record = [log_val_record; new_log_val];
%     
%     if isnan(new_log_val)
%         pause
%     end
        
        if norm(new_z_est(:) - z_est(:)) <= 1e-4
            outer_n_iter = kk;
            break;
        else
            h_est = new_h_est;
            lambda_est = new_lambda_est;
            pi_est = new_pi_est;
            z_est = new_z_est;
%             log_val = new_log_val;
        end
        
%         sigma_est = sqrt(1./lambda_est);

end

%% record the converged vals
    record.h_est = new_h_est;
    record.lambda_est = new_lambda_est;
    record.z_est = new_z_est;
    record.pi_est = new_pi_est;
    record.inner_n_iter_record = inner_n_iter_record;
%     record.log_val = new_log_val;
%     record.log_val_record = log_val_record;
    if outer_n_iter~=0
        record.outer_n_iter = outer_n_iter;
    else
        record.outer_n_iter = n_iter;
    end
    
    record.max_n_pt = max_n_pt;
    record.max_n_s = max_n_s;
    
    record.log_alpha = log_alpha;
    record.delta = inner_record.delta;
    record.zeta = inner_record.zeta;
    record.rho = inner_record.rho;
    record.eta = inner_record.eta;
    record.tau = inner_record.tau;
    