function z_est_online = func_ci_mul_pre_bias_online(X, pi_est, h_est, lambda_est)

% num of states
[~, n_s] = size(X);
% num of diff levels
n_l = size(h_est,2);
z_est_online = zeros(1,n_s);

for j = 1:n_s
   idx_j = ~isnan(X(:,j));
   cur_mu = median(X(idx_j, j));
   cur_nu = 0.5*1/(var(X(idx_j, j)) + 1e-4);
   
   cur_n_pt = sum(idx_j);
   cur_h = zeros(cur_n_pt+1, n_l);
   cur_lambda = zeros(cur_n_pt+1, n_l);
   
   for k = 1:n_l
       cur_h(:,k) = [cur_mu; X(idx_j,j)-h_est(idx_j,k)]; % col vec
       cur_lambda(:,k) = [cur_nu; lambda_est(idx_j,k)]; % col vec
       z_est_online(j) = z_est_online(j) + pi_est(k)*(cur_h(:,k).'*cur_lambda(:,k))/sum(cur_lambda(:,k));    
   end
   
end
