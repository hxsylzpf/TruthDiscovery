function record_mpm_mul = func_em_ci_mul_pre_bias_ini_mul_run(n_pt, n_s, n_l, X, ini_z, ini_lambda, ini_h, prior_mu, prior_nu, ...
    prior_mu_h, prior_nu_h, prior_a_mul, prior_b_mul)

n_run = 1; % 100
record_max_log_val = zeros(n_run,1);
record_z_val = zeros(n_run, n_s);
cpu_time = zeros(n_run,1);

for i = 1:n_run
    if i == 1
        cur_ini_z = ini_z;
        % ini lambda will not be used
        cur_ini_lambda = ini_lambda;
    else
%         rng(i);
        % pertube ini_z
        cur_ini_z = ini_z.*(0.5 + rand(1,n_s));
        % ini lambda will not be used
        cur_ini_lambda = ini_lambda.*(0.5 + rand(n_pt, n_l));
    end
    
    tic
    record_mul_run(i) = func_em_ci_mul_pre_bias_ini(n_pt, n_s, n_l, X, cur_ini_z, cur_ini_lambda, ini_h, prior_mu, prior_nu, ...
        prior_mu_h, prior_nu_h, prior_a_mul, prior_b_mul);
    cpu_time(i) = toc;
    record_max_log_val(i) = record_mul_run(i).log_val;
    record_z_val(i,:) = record_mul_run(i).z_est;
end

[max_val, idx] = max(record_max_log_val);

% max_val

% these first two records are for the best run
record_mpm_mul = record_mul_run(idx);
record_mpm_mul.cpu_time = mean(cpu_time);

% these are for the debug purpose, which records the results in all the runs
record_mpm_mul.record_max_log_val = record_max_log_val;
record_mpm_mul.record_z_val = record_z_val;
