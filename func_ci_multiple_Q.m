function log_val = func_ci_multiple_Q(n_pt, n_s, n_l, X, pi, r, lambda, z, prior_mu, prior_nu, prior_a, prior_b)

log_val = 0;

for j = 1:n_s
    log_val = log_val + 0.5*mylog(prior_nu(j)) - 0.5*prior_nu(j)*(z(j) - prior_mu(j))^2;
    for k = 1:n_l
        temp = mylog(pi(k));
        for i = 1:n_pt
            if ~isnan(X(i,j))
                temp = temp + 0.5*mylog(lambda(i,k)) - 0.5*lambda(i,k)*(X(i,j) - z(j))^2;
            end
        end
        log_val = log_val + r(j,k)*temp;
    end
end

for i = 1:n_pt
    for k = 1:n_l
        log_val = log_val + (prior_a(i,k)-1)*mylog(lambda(i,k)) - prior_b(i,k)*lambda(i,k);
    end
end
