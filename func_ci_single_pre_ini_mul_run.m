function record_tbp_mul = func_ci_single_pre_ini_mul_run(n_pt, n_s, X, ini_z, prior_mu, prior_nu, prior_a, prior_b)

% no need to run 100 times, just once is enough
n_run = 1; % 100
record_max_log_val = zeros(n_run,1);
record_z_val = zeros(n_run, n_s);
cpu_time = zeros(n_run,1);

for i = 1:n_run
    if i==1
        cur_ini_z = ini_z;
    else
%         rng(i);
        cur_ini_z = ini_z.*(0.5 + rand(1,n_s));
    end
    
%     cur_ini_z
    
    tic
    record_mul_run(i) = func_ci_single_pre_ini...
        (n_pt, n_s, X, cur_ini_z, prior_mu, prior_nu, prior_a, prior_b);
    cpu_time(i) = toc;
    record_max_log_val(i) = record_mul_run(i).log_val;
    record_z_val(i,:) = record_mul_run(i).z_est;
end

[max_val, idx] = max(record_max_log_val);

max_val

% these first two records are for the best run
record_tbp_mul = record_mul_run(idx);
record_tbp_mul.cpu_time = mean(cpu_time);

% these are for the debug purpose, which records the results in all the runs
record_tbp_mul.record_max_log_val = record_max_log_val;
record_tbp_mul.record_z_val = record_z_val;
