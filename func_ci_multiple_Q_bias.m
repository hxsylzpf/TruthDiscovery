function log_val = func_ci_multiple_Q_bias(n_pt, n_s, n_l, X, pi, r, h, lambda, z, ...
    prior_mu, prior_nu, prior_mu_h, prior_nu_h, prior_a, prior_b)

log_val = 0;

for j = 1:n_s
    log_val = log_val + 0.5*mylog(prior_nu(j)) - 0.5*prior_nu(j)*(z(j) - prior_mu(j))^2;
    for k = 1:n_l
        % initialization
        temp = log(pi(k));
        for i = 1:n_pt
            if ~isnan(X(i,j))
                temp = temp + 0.5*log(lambda(i,k)) - 0.5*lambda(i,k)*(X(i,j) - z(j) - h(i,k))^2;
            end
        end
        log_val = log_val + r(j,k)*temp;
    end
end

for i = 1:n_pt
    for k = 1:n_l
        % prior on h
%         log_val = log_val + 0.5*log(prior_nu_h(i,k)) - 0.5*prior_nu_h(i,k)*(h(i,k) - prior_mu_h(i,k))^2;
        % prior on lambda
        log_val = log_val + (prior_a(i,k)-1)*log(lambda(i,k)) - prior_b(i,k)*lambda(i,k);
    end
end
