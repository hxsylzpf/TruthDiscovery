function [pi, mu, sigma] = incre_mix_gau(data, n_l, pi_est, mu_est, sigma_est)

n_data = size(data,1);

s1 = zeros(1, n_l);
s2 = zeros(1, n_l);
s3 = zeros(1, n_l);
    
for i = 1:n_data
%% E-step
    gamma = zeros(1, n_l);
    
    step = i^(-0.7);
    
    for k = 1:n_l
        gamma(k) = pi_est(k)*my_normpdf(data(i),mu_est(k),sigma_est(k));
    end

    norm_gamma = gamma/sum(gamma);
    
    for k = 1:n_l
        s1(k) = (1-step)*s1(k) + step*norm_gamma(k);
        s2(k) = (1-step)*s2(k) + step*norm_gamma(k)*data(i);
        s3(k) = (1-step)*s3(k) + step*norm_gamma(k)*data(i)^2;
    end
    
%% M-step
if i > 30
    % directly modify
    for k = 1:n_l
        pi_est(k) = s1(k);
        mu_est(k) = s2(k)/s1(k);
        sigma_est(k) = sqrt((s3(k) - 2*s2(k)*mu_est(k) + s1(k)*mu_est(k)^2)/s1(k));
    end
end

end

pi = pi_est;
mu = mu_est;
sigma = sigma_est;
