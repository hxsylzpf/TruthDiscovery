function z_est_online = func_ci_single_pre_bias_online(x, h_est, lambda_est)

% x is for a single event
% h_est and lambda_est are also selected out

    cur_mu = median(x);
    cur_nu = 1e-6*1/(var(x - cur_mu) + 1e-4);
    cur_h = [cur_mu; x-h_est]; % col vec
    cur_lambda = [cur_nu; lambda_est]; % col vec
    z_est_online = (cur_h.'*cur_lambda)/sum(cur_lambda);
