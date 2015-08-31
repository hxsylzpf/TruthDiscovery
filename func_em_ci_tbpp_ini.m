%% model parameter estimation from given data - EM
% tbpp - truth, bias, precision and preference
function record = func_em_ci_tbpp_ini(n_pt, n_s, n_l, X, ini_z, n_clu_per_pt, ini_pi, ini_h, ini_lambda, prior_mu, prior_nu, prior_a, prior_b, style)
% n_pt - num of pts
% n_s - num of states/targets
% n_l - num of levels
% parameter initialization is very important for EM
% bad initialization % can result in bad estimation

%% initialization
n_iter = 100;
outer_n_iter = 0;

log_val = -inf;
log_val_record =[];
inner_n_iter_record = [];


%% if the input ini_pi, ini_h and ini_lambda are in the matrix format
if strcmp(style, 'matrix')

% initialize para values
z_est = ini_z;
pi_est = ini_pi;
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
                for k = 1:n_l
                    temp_zeta(i,j,k) = pi_est(i,k)*my_normpdf(X(i,j),z_est(j)+h_est(i,k),sqrt(1/lambda_est(i,k)));
                end
                zeta(i,j,:) = temp_zeta(i,j,:)/sum(temp_zeta(i,j,:));
            end
        end
    end
        
    %% M-step - estimate model paras
    new_pi_est = zeros(n_pt, n_l);
    for i = 1:n_pt
        idx_i = ~isnan(X(i,:));
        Ni = sum(idx_i);
        for k = 1:n_l
            new_pi_est(i,k) = (sum(zeta(i,idx_i,k))+0.001)/(Ni+0.001*n_l);
        end
    end    
        
    %% the last updated z value is used as the initial values for para reestimation
    inner_record = func_solve_z_h_lambda_tbpp(n_pt, n_s, n_l, X, zeta, z_est, prior_mu, prior_nu, prior_a, prior_b);
    
    new_z_est = inner_record.z_est;
    new_h_est = inner_record.h_est;
    new_lambda_est = inner_record.lambda_est;
        
    new_log_val = func_tbpp_log_likelihood(n_pt, n_s, n_l, X, pi_est, h_est, lambda_est, z_est, ...
    prior_mu, prior_nu, prior_a, prior_b, n_clu_per_pt, style);
    
    log_val_record = [log_val_record; new_log_val];
    
        if new_log_val - log_val < 1e-4
            outer_n_iter = kk;
            break;
        else
            h_est = new_h_est;
            lambda_est = new_lambda_est;
            pi_est = new_pi_est;
            z_est = new_z_est;
            log_val = new_log_val;
        end
        
end

%% if the input ini_pi, ini_h and ini_lambda are in the cell format
elseif strcmp(style, 'cell')

    
% initialize para values
z_est = ini_z;
pi_est = ini_pi;
h_est = ini_h;
lambda_est = ini_lambda;

for kk = 1:n_iter
    
    %% E-step - calculate expectation
    %% zeta
    zeta = cell(n_pt, n_s);
    for i = 1:n_pt
        for j = 1:n_s
            if ~isnan(X(i,j))
                if n_clu_per_pt(i) == 1
                    zeta{i,j} = 1;
                else
                    temp_zeta = zeros(1, n_clu_per_pt(i));
                    for k = 1:n_clu_per_pt(i)
                        temp_zeta(k) = pi_est{i}(k)*my_normpdf(X(i,j),z_est(j)+h_est{i}(k),sqrt(1/lambda_est{i}(k)));
                    end
                end
                zeta{i,j} = temp_zeta/sum(temp_zeta);
            end
        end
    end
        
    %% M-step - estimate model paras
    new_pi_est = cell(n_pt, 1);
    for i = 1:n_pt
        idx_i = ~isnan(X(i,:));
        Ni = sum(idx_i);
        if n_clu_per_pt(i) == 1
            new_pi_est{i} = 1;
        else
            for k = 1:n_clu_per_pt(i)
                temp = 0;
                for j = 1:n_s
                    if ~isnan(X(i,j))
                        temp = temp + zeta{i,j}(k);
                    end
                end
                new_pi_est{i}(k) = temp/Ni;
            end
        end
    end    
        
    %% the last updated z value is used as the initial values for para reestimation
    inner_record = func_solve_z_h_lambda_tbpp_cell(n_pt, n_s, X, n_clu_per_pt, zeta, z_est, prior_mu, prior_nu, prior_a, prior_b);
    
    new_z_est = inner_record.z_est;
    new_h_est = inner_record.h_est;
    new_lambda_est = inner_record.lambda_est;
        
    new_log_val = func_tbpp_log_likelihood(n_pt, n_s, n_l, X, pi_est, h_est, lambda_est, z_est, ...
    prior_mu, prior_nu, prior_a, prior_b, n_clu_per_pt, style);
    
    log_val_record = [log_val_record; new_log_val];
    
        if new_log_val - log_val < 1e-4
            outer_n_iter = kk;
            break;
        else
            h_est = new_h_est;
            lambda_est = new_lambda_est;
            pi_est = new_pi_est;
            z_est = new_z_est;
            log_val = new_log_val;
        end
        
end   
    



end

    record.h_est = new_h_est;
    record.lambda_est = new_lambda_est;
    record.z_est = new_z_est;
    record.pi_est = new_pi_est;
    record.inner_n_iter_record = inner_n_iter_record;
    record.log_val = new_log_val;
    record.log_val_record = log_val_record;
    if outer_n_iter~=0
        record.outer_n_iter = outer_n_iter;
    else
        record.outer_n_iter = n_iter;
    end
    