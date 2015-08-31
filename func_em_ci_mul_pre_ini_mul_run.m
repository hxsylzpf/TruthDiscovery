function record_mpm_mul = func_em_ci_mul_pre_ini_mul_run(n_pt, n_s, n_l, X, prior_mu, prior_nu, prior_a_mul, prior_b_mul)

n_run = 50;
record_log_val = zeros(n_run,1);

for i = 1:n_run
    record_mul_run(i) = func_em_ci_mul_pre_ini(n_pt, n_s, n_l, X, prior_mu, prior_nu, prior_a_mul, prior_b_mul);
    record_log_val(i) = record_mul_run(i).log_val;
end

[max_val, idx] = max(record_log_val);

max_val

record_mpm_mul = record_mul_run(idx);