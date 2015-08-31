function log_val = func_tbpp_log_likelihood(n_pt, n_s, n_l, X, pi_est, h_est, lambda_est, z_est, ...
    prior_mu, prior_nu, prior_a, prior_b, n_clu_per_pt, style)

log_val = 0;

if strcmp(style, 'matrix')
   for j = 1:n_s
       % prior
        log_val = log_val + 0.5*mylog(prior_nu(j)) - 0.5*prior_nu(j)*(z_est(j) - prior_mu(j))^2;
        for i = 1:n_pt
            if ~isnan(X(i,j))
                temp = 0;
                for k = 1:n_l
                    temp = temp + pi_est(i,k)*my_normpdf(X(i,j),z_est(j)+h_est(i,k),sqrt(1/lambda_est(i,k)));
                end
                log_val = log_val + log(temp);
            end
        end
    end

    for i = 1:n_pt
        for k = 1:n_l
            % prior on lambda
            log_val = log_val + (prior_a(i,k)-1)*log(lambda_est(i,k)) - prior_b(i,k)*lambda_est(i,k);
        end
    end

    
elseif strcmp(style, 'cell')
    
   for j = 1:n_s
       % prior
        log_val = log_val + 0.5*mylog(prior_nu(j)) - 0.5*prior_nu(j)*(z_est(j) - prior_mu(j))^2;
        for i = 1:n_pt
            if ~isnan(X(i,j))
                temp = 0;
                for k = 1:n_clu_per_pt(i)
                    temp = temp + pi_est{i}(k)*my_normpdf(X(i,j),z_est(j)+h_est{i}(k),sqrt(1/lambda_est{i}(k)));
                end
                log_val = log_val + log(temp);
            end
        end
    end

    for i = 1:n_pt
        for k = 1:n_clu_per_pt(i)
            % prior on lambda
            log_val = log_val + (prior_a{i}(k)-1)*log(lambda_est{i}(k)) - prior_b{i}(k)*lambda_est{i}(k);
        end
    end
    
    
end