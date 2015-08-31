function record = func_ci_single_pre_bias_online(X, z, scale_fac, h_est, lambda_est, n_est_range)
% n_test is the num of times you test your online estimator
% n_est_range is the num of estimations per target
% (note: n_pt >= n_est as some pt may not participate in a particular task)

N = length(n_est_range);
n_s = size(X,2);
err_record = zeros(N,n_s);
rmse_record = zeros(N,1);
cpu_time_record = zeros(N,1);

for k = 1:N
    n_est = n_est_range(k);
    z_est_online = zeros(1,n_s);
    time_online = zeros(n_s,1);
    % test on all the events
    for j = 1:n_s
       % find who participated in this task
       idx_j = ~isnan(X(:,j));
       % randomly pick out n_est pts
       out_idx = rand_select_idx(idx_j, 1, n_est);
       tic
       cur_mu = median(X(out_idx, j));
       cur_nu = 0.01*1/(var(X(out_idx, j)) + 1e-4);
       cur_h = [cur_mu; X(out_idx,j)-h_est(out_idx)]; % col vec
       cur_lambda = [cur_nu; lambda_est(out_idx)]; % col vec
       z_est_online(j) = (cur_h.'*cur_lambda)/sum(cur_lambda);
       time_online(j) = toc;
    end
    err_record(k,:) = round(z_est_online/scale_fac)-(z/scale_fac).';
    rmse_record(k) = my_rmse(err_record(k,:));
    cpu_time_record(k) = mean(time_online);
end

record.err_record = err_record;
record.rmse_record = rmse_record;
record.cpu_time_record = cpu_time_record;
