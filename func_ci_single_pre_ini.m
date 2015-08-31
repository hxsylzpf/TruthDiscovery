%% model parameter estimation from given data
function record = func_ci_single_pre_ini(n_pt, n_s, X, ini_z, prior_mu, prior_nu, prior_a, prior_b)
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% initialization
n_iter = 100;
% lambda_est = ones(n_pt, 1);
z_est = ini_z;
log_val = -inf;
log_val_record = [];

for kk = 1:n_iter

    new_lambda_est = zeros(n_pt, 1);
    %% estimate lambda
    for i = 1:n_pt
        idx_i = (~isnan(X(i,:)));
        ni = sum(idx_i);
        
%         % print out to check for errors
%         v1 = X(i,idx_i) - z_est(idx_i)
%         v2 = sum((X(i,idx_i) - z_est(idx_i)).^2)

        new_lambda_est(i) = (ni/2 + prior_a(i) - 1)/(1/2*sum((X(i,idx_i) - z_est(idx_i)).^2) + prior_b(i));
    end
    
%     new_lambda_est
    
    new_log_val = func_ci_single_Q(n_pt, n_s, X, new_lambda_est, z_est, prior_mu, prior_nu, prior_a, prior_b)
    
    log_val_record = [log_val_record; new_log_val];
    
        if new_log_val - log_val < 1e-5
            break;
        else
            lambda_est = new_lambda_est;
            log_val = new_log_val;
        end

%         sigma_est = sqrt(1./lambda_est);
        
    %% estimate z
    for j = 1:n_s
        idx_j = (~isnan(X(:,j)));
        z_est(j) = (prior_nu(j)*prior_mu(j) + lambda_est(idx_j).'*X(idx_j,j))/(prior_nu(j)+sum(lambda_est(idx_j)));
%         z_est(j) = (lambda_est(idx_j).'*X(idx_j,j))/(sum(lambda_est(idx_j)));
    end  
    
    
end
   
    % after convergence
    data_ll = func_ci_single_data_ll(n_pt, n_s, X, lambda_est, z_est);

    record.lambda_est = new_lambda_est;
%     record.z_est = round(z_est);
    record.z_est = z_est;
    record.data_ll = data_ll;    
    record.log_val = new_log_val;
    record.log_val_record = log_val_record;
    