%% model parameter estimation from given data
function record = func_ci_single_pre_bias_ini(n_pt, n_s, X, ini_z, ini_h, prior_mu, prior_nu, prior_mu_h, prior_nu_h, prior_a, prior_b)
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% initialization
n_iter = 100;
z_est = ini_z;
h_est = ini_h;
% h_est = ini_h;
% h_est = zeros(n_pt, 1); % -0.5 + rand
log_val = -inf;
log_val_record = [];

for kk = 1:n_iter

    new_lambda_est = zeros(n_pt, 1);
    new_h_est = zeros(n_pt, 1);
    new_z_est = zeros(size(z_est));
    
    %% estimate h and lambda
    for i = 1:n_pt
        idx_i = ~isnan(X(i,:));
        ni = sum(idx_i);
%         new_h_est(i) = (prior_mu_h(i)*prior_nu_h(i) + new_lambda_est(i)*sum(X(i, idx_i) - z_est(idx_i)))/...
%             (new_lambda_est(i)*ni + prior_nu_h(i));
        new_h_est(i) = (sum(X(i, idx_i) - z_est(idx_i)))/ni;
    end
    
    for i = 1:n_pt
        idx_i = ~isnan(X(i,:));
        ni = sum(idx_i);
        new_lambda_est(i) = (ni/2 + prior_a(i) - 1)/(1/2*sum((X(i,idx_i) - z_est(idx_i) - new_h_est(i)).^2) + prior_b(i));
%         new_lambda_est(i) = ni/sum((X(i,idx_i) - z_est(idx_i) - new_h_est(i)).^2);
    end
        

    
    %% estimate z
    for j = 1:n_s
        idx_j = ~isnan(X(:,j));
        new_z_est(j) = (prior_nu(j)*prior_mu(j) + new_lambda_est(idx_j).'*(X(idx_j,j) - new_h_est(idx_j)))/(prior_nu(j)+sum(new_lambda_est(idx_j)));
%         new_z_est(j) = (new_lambda_est(idx_j).'*(X(idx_j,j) - new_h_est(idx_j)))/sum(new_lambda_est(idx_j));
    end
    
    %% evaluate obj func value
    new_log_val = func_ci_single_Q_bias(n_pt, n_s, X, new_h_est, new_lambda_est, new_z_est, prior_mu, prior_nu, prior_mu_h, prior_nu_h, prior_a, prior_b);
    log_val_record = [log_val_record; new_log_val];
    
    if new_log_val - log_val < 1e-5 % norm(new_lambda_est(:) - lambda_est(:))
       break;
    else
%         h_est = new_h_est;
        h_est = new_h_est;
        z_est = new_z_est;
        log_val = new_log_val;
    end
        
%         if norm(new_lambda_est - lambda_est) + norm(new_h_est - h_est)< 1e-4
%             break;
%         else
%             lambda_est = new_lambda_est;
%             h_est = new_h_est;
%             z_est = new_z_est;
%         end

%         sigma_est = sqrt(1./lambda_est);

end

    % after convergence
    data_ll = func_ci_single_bias_data_ll(n_pt, n_s, X, new_h_est, new_lambda_est, new_z_est);

    record.h_est = new_h_est;
    record.lambda_est = new_lambda_est;
    record.z_est = new_z_est;
    record.data_ll = data_ll;
    record.log_val = new_log_val;
    record.log_val_record = log_val_record;
    