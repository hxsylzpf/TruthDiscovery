function [pi, mu, sigma] = batch_mix_gau(data, n_l, pi_est, mu_est, sigma_est)

n_data = size(data,1);

gamma = zeros(n_data, n_l);

for iter = 1:50
%% E-step
for i = 1:n_data
    for k = 1:n_l
        gamma(i,k) = pi_est(k)*my_normpdf(data(i),mu_est(k),sigma_est(k));        
    end
end

norm_gamma = prob_mat_nlz(gamma, 'row');

%% M-step

N = zeros(1,n_l);
new_pi_est = zeros(size(pi_est));
new_mu_est = zeros(size(mu_est));
new_sigma_est = zeros(size(sigma_est));

for k = 1:n_l
    N(k) = sum(norm_gamma(:,k));
    new_pi_est(k) = N(k)/n_data;
    new_mu_est(k) = sum(norm_gamma(:,k).*data)/N(k);
    new_sigma_est(k) = sqrt(sum(norm_gamma(:,k).*(data - mu_est(k)).^2)/N(k));
end

if norm(new_pi_est - pi_est) + norm(new_mu_est - mu_est) + norm(new_sigma_est - sigma_est) < 1e-3
    break;
else
    pi_est = new_pi_est;
    mu_est = new_mu_est;
    sigma_est = new_sigma_est;
end

end

pi = pi_est;
mu = mu_est;
sigma = sigma_est;
