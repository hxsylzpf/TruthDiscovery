%% model parameter estimation from given data - EM
% tbpp with hidden diff - truth, bias, precision and perception
function record = func_em_ci_tbpp_hidden_ini(n_pt, n_s, n_l, X, per_diff, ini_z, ini_eta, ini_w, ini_w0, ini_h, ini_lambda,  ...
    prior_mu, prior_nu, prior_a, prior_b)
% n_pt - num of pts
% n_s - num of states/targets
% n_l - num of levels
% parameter initialization is very important for EM
% bad initialization % can result in bad estimation

%% initialization
n_iter = 50;
outer_n_iter = 0;

% log_val = -inf;
% log_val_record =[];
inner_n_iter_record = [];

% initialize para values
z_est = ini_z;
eta_est = ini_eta;
w_est = ini_w;
w0_est = ini_w0;
h_est = ini_h;
lambda_est = ini_lambda;

for kk = 1:n_iter
    
    %% E-step - calculate expectation
    %% zeta
    temp_zeta = zeros(n_pt, n_s, n_l);
    zeta = zeros(n_pt, n_s, n_l);
    for i = 1:n_pt
        for j = 1:n_s
            if ~isnan(X(i,j))
                [~, norm_temp_est] = func_logistic_reg_mul_pre(n_l, w_est(i,:), w0_est(i,:), z_est(j));
                for k = 1:n_l
                    cur_I = (per_diff(i,j) == k);
                    temp_zeta(i,j,k) = eta_est(i)^cur_I*(1-eta_est(i))^(1-cur_I) ...
                        *my_normpdf(X(i,j),z_est(j)+h_est(i,k),sqrt(1/lambda_est(i,k))) ...
                        *norm_temp_est(k);
                end
                zeta(i,j,:) = temp_zeta(i,j,:)/sum(temp_zeta(i,j,:));
            end
        end
    end
        
    %% M-step - estimate model paras
    %% eta
    new_eta_est = zeros(size(eta_est));
    for i = 1:n_pt
        idx = ~isnan(X(i,:));
        ni = sum(idx);
        for j = 1:n_s
            if ~isnan(X(i,j))
                new_eta_est(i) = new_eta_est(i) + zeta(i,j,per_diff(i,j));
            end
        end
        new_eta_est(i) = new_eta_est(i)/ni;
    end

    %% w, w0
    new_w_est = zeros(size(w_est));
    new_w0_est = zeros(size(w0_est));
    
    for i = 1:n_pt
        idx = ~isnan(per_diff(i,:));
        cur_z = z_est(idx);
        cur_zeta = zeta(i,idx,:);

        [cur_w_est, cur_w0_est] = func_tbpp_hidden_w(n_l, cur_zeta, cur_z.');
        new_w_est(i,:) = cur_w_est.';
        new_w0_est(i,:) = cur_w0_est.';        
    end
    
    %% the last updated z value is used as the initial values for para reestimation
    inner_record = func_solve_z_h_lambda_tbpp(n_pt, n_s, n_l, X, zeta, z_est, prior_mu, prior_nu, prior_a, prior_b);
    
    new_z_est = inner_record.z_est;
    new_h_est = inner_record.h_est;
    new_lambda_est = inner_record.lambda_est;
        
%     new_log_val = func_tbpp_log_likelihood(n_pt, n_s, n_l, X, pi_est, h_est, lambda_est, z_est, ...
%     prior_mu, prior_nu, prior_a, prior_b, n_clu_per_pt, style);
%     
%     log_val_record = [log_val_record; new_log_val];
    
        if norm(new_z_est - z_est) < 1e-3 % new_log_val - log_val < 1e-4
            outer_n_iter = kk;
            break;
        else
            eta_est = new_eta_est;
            w_est = new_w_est;
            w0_est = new_w0_est;
            h_est = new_h_est;
            lambda_est = new_lambda_est;
            z_est = new_z_est;
%             log_val = new_log_val;
        end
          
end

    record.eta_est = new_eta_est;
    record.w_est = new_w_est;
    record.w0_est = new_w0_est;
    record.h_est = new_h_est;
    record.lambda_est = new_lambda_est;
    record.z_est = new_z_est;
    record.inner_n_iter_record = inner_n_iter_record;
%     record.log_val = new_log_val;
%     record.log_val_record = log_val_record;
    if outer_n_iter~=0
        record.outer_n_iter = outer_n_iter;
    else
        record.outer_n_iter = n_iter;
    end
    