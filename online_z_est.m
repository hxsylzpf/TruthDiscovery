function z_est = online_z_est(n_l, x, lambda_est, prior_mean, prior_pre, pi_est)

z_est = 0;

idx_i = ~isnan(x);

for k = 1:n_l
    temp1 = [prior_mean x(idx_i)];
    temp2 = [prior_pre;lambda_est(idx_i,k)];
    z_est = z_est + pi_est(k)*(temp1*temp2)/sum(temp2);
    
end
