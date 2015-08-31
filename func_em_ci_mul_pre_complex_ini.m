%% model parameter estimation from given data - EM
function record = func_em_ci_mul_pre_complex_ini(n_pt, n_s, n_l, X, prior_mu, prior_nu, prior_a, prior_b)
% n_s - num of states
% n_l - num of levels
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% initialization
n_iter = 100;
pi_est = 1/n_l*ones(n_l, 1);
lambda_est = rand(n_pt, n_l);

for kk = 1:n_iter

    %% E-step - calculate expectation
    %% alpha
    % use log to calculate to avoid num underflow
    temp_log_alpha = zeros(n_s, n_l);
    for j = 1:n_s
        for k = 1:n_l
            temp_log_alpha(j,k) = log(pi_est(k));
            for i = 1:n_pt
               if ~isnan(X(i,j))
                   temp_log_alpha(j,k) = temp_log_alpha(j,k) + log(my_normpdf(X(i,j),prior_mu(j),sqrt(1/prior_nu(j)+1/lambda_est(i,k))));
               end
            end
        end
    end
    
    % convert back to num
    temp_alpha = exp(temp_log_alpha);
    % normalize
    alpha = prob_mat_nlz(temp_alpha, 'row');

    %% w1, w2
    w1 = zeros(n_s, n_l);
    w2 = zeros(n_s, n_l);
    for j = 1:n_s
        idx_j = ~isnan(X(:,j));
        for k = 1:n_l
            w1(j,k) = (prior_nu(j)*prior_mu(j) + lambda_est(idx_j,k).'*X(idx_j,j))/(prior_nu(j)+sum(lambda_est(idx_j,k)));
            w2(j,k) = 1/(prior_nu(j)+sum(lambda_est(idx_j,k))) + w1(j,k)^2;
        end
    end

    %% beta, gamma
    beta = alpha.*w1;
    gamma = alpha.*w2;
    
    %% M-step - estimate model paras
    %% pi
    new_pi_est = sum(alpha)/n_s;
    new_pi_est = new_pi_est.';
    new_lambda_est = zeros(size(lambda_est));
    
    %% lambda
    for i = 1:n_pt
        idx_i = ~isnan(X(i,:));
        for k = 1:n_l
            new_lambda_est(i,k) = (1/2*sum(alpha(idx_i, k)) + prior_a(i,k) - 1)/...
                (1/2*(X(i,idx_i).^2*alpha(idx_i, k) - 2*X(i,idx_i)*beta(idx_i,k) + sum(gamma(idx_i,k))) + prior_b(i,k));
        end
    end
    
        if norm(new_lambda_est(:) - lambda_est(:)) + norm(new_pi_est - pi_est)< 1e-4
            break;
        else
            lambda_est = new_lambda_est;
            pi_est = new_pi_est;
        end
        
        sigma_est = sqrt(1./lambda_est);

end

    z_est = sum(beta, 2).';

    record.lambda_est = lambda_est;
%     record.z_est = round(z_est);
    record.z_est = z_est;
    
    