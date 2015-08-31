%% model parameter estimation from given data - EM
function record = func_tbp_ori(n_pt, n_s, n_l, X, ini_z, ini_h, ini_lambda, prior_mu, prior_nu, prior_a, prior_b)
% n_s - num of states
% n_l - num of levels
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% initialization
n_iter = 100;
pi_est = 1/n_l*ones(1, n_l);

log_val = -inf;
% initialize para values
z_est = ini_z;
h_est = ini_h;
lambda_est = ini_lambda;

for kk = 1:n_iter

%     z_est
    
    %% E-step - calculate expectation
    %% alpha
    alpha = zeros(n_s, n_l);
    
    temp1 = zeros(n_s, n_l);
    temp2 = zeros(n_s, n_l);
    
   %% for all targets
    for j = 1:n_s
        idx_j = ~isnan(X(:,j));
        
        % for all level
        for k = 1:n_l            
            temp1(j,k) = lambda_est(idx_j,k).'*(X(idx_j,j) - h_est(idx_j,k));
            temp2(j,k) = sum(lambda_est(idx_j,k));
            
            temp_mu = temp1(j,k)/temp2(j,k);
            temp_nu = temp2(j,k);
            
            alpha(j,k) = pi_est(k)*normpdf(temp_mu,z_est(j),sqrt(1/temp_nu));
        end    
        
        alpha(j,:) = alpha(j,:)/sum(alpha(j,:));
    end
%     pi_est
%     alpha
%     sum(alpha)

    %% M-step - estimate model paras
    
    new_pi_est = (sum(alpha)+1)/(n_s+n_l);
    
    %% the last updated z value is used as the initial values for para reestimation
    inner_record = func_solve_z_h_lambda(n_pt, n_s, n_l, X, alpha, z_est, lambda_est, prior_mu, prior_nu, [], [], prior_a, prior_b);
    
    new_z_est = inner_record.z_est;
    new_h_est = inner_record.h_est;
    new_lambda_est = inner_record.lambda_est;
    
    new_log_val = func_ci_multiple_Q_bias(n_pt, n_s, n_l, X, new_pi_est, alpha, new_h_est, new_lambda_est, new_z_est, ...
    prior_mu, prior_nu, [], [], prior_a, prior_b);
    
        if  new_log_val - log_val < 1e-4 % norm(new_h_est(:) - h_est(:)) + norm(new_lambda_est(:) - lambda_est(:)) < 1e-3
            break;
        else
            h_est = new_h_est;
            lambda_est = new_lambda_est;
            pi_est = new_pi_est;
            z_est = new_z_est;
            log_val = new_log_val;
        end
        

end

    % after convergence
    
    record.h_est = new_h_est;
    record.lambda_est = new_lambda_est;
    record.z_est = new_z_est;
    record.pi_est = new_pi_est;
    