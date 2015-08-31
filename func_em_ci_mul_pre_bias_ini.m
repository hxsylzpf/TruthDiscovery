%% model parameter estimation from given data - EM
function record = func_em_ci_mul_pre_bias_ini(n_pt, n_s, n_l, X, ini_z, ini_lambda, ini_h, prior_mu, prior_nu, ....
    prior_mu_h, prior_nu_h, prior_a, prior_b)
% n_s - num of states
% n_l - num of levels
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% initialization
n_iter = 100;
pi_est = 1/n_l*ones(1, n_l);
outer_n_iter = 0;

log_val = -inf;
log_val_record =[];
inner_n_iter_record = [];

% initialize para values
z_est = ini_z;
% [z_est, ~, ~] = ini_z_lambda(X);
lambda_est = ini_lambda; % rand(n_pt, n_l);
h_est = ini_h; % 0.1*(-0.5+rand(n_pt, n_l));

for kk = 1:n_iter

%     z_est
    
    %% E-step - calculate expectation
    %% alpha
    % use log to calculate to avoid num underflow
    temp_log_alpha = zeros(n_s, n_l);
    for j = 1:n_s
        for k = 1:n_l
            temp_log_alpha(j,k) = log(pi_est(k));
            for i = 1:n_pt
               if ~isnan(X(i,j))
                   temp_log_alpha(j,k) = temp_log_alpha(j,k) + log(my_normpdf(X(i,j),z_est(j)+h_est(i,k),sqrt(1/lambda_est(i,k))));
               end
            end
        end
    end
        
    % convert back to num, also prevent numerical error
    temp_alpha = exp(temp_log_alpha);
    % normalize
    alpha = prob_mat_nlz(temp_alpha, 'row');
    
%     pi_est
%     alpha
%     sum(alpha)

    %% M-step - estimate model paras
    
    new_pi_est = (sum(alpha)+1)/(n_s+n_l);
    
    %% the last updated z value is used as the initial values for para reestimation
    inner_record = func_solve_z_h_lambda(n_pt, n_s, n_l, X, alpha, z_est, lambda_est, prior_mu, prior_nu, prior_mu_h, prior_nu_h, prior_a, prior_b);
    
    new_z_est = inner_record.z_est;
    new_h_est = inner_record.h_est;
    new_lambda_est = inner_record.lambda_est;
%     inner_n_iter_record = [inner_n_iter_record record.n_iter_to_converge];
   
    new_log_val = func_ci_multiple_Q_bias(n_pt, n_s, n_l, X, new_pi_est, alpha, new_h_est, new_lambda_est, new_z_est, ...
        prior_mu, prior_nu, prior_mu_h, prior_nu_h, prior_a, prior_b);

    log_val_record = [log_val_record; new_log_val];
    
%     if isnan(new_log_val)
%         alpha
%         pause
%     end
        
        if  new_log_val - log_val < 1e-4 % norm(new_lambda_est(:) - lambda_est(:))
            outer_n_iter = kk;
            break;
        else
            h_est = new_h_est;
            lambda_est = new_lambda_est;
            pi_est = new_pi_est;
            z_est = new_z_est;
            log_val = new_log_val;
        end
        
%         sigma_est = sqrt(1./lambda_est);

end

    % after convergence
    data_ll = func_ci_multiple_bias_incmpl_logll(n_pt, n_s, n_l, X, new_pi_est, new_h_est, new_lambda_est, new_z_est);
    
    record.h_est = new_h_est;
    record.lambda_est = new_lambda_est;
    record.z_est = new_z_est;
    record.pi_est = new_pi_est;
    [~, record.r_est] = max(alpha, [], 2);
    record.inner_n_iter_record = inner_n_iter_record;
    record.log_val = new_log_val;
    record.log_val_record = log_val_record;
    record.data_ll = data_ll;
    
    if outer_n_iter~=0
        record.outer_n_iter = outer_n_iter;
    else
        record.outer_n_iter = n_iter;
    end
    