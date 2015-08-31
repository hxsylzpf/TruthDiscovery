function val = func_gmm_pdf(x, K, pi, mu, sigma)

% x is the data
% K is the num of components
% sigma is the std, not variance

val = 0;
for i = 1:K
    val = val + pi(i)*normpdf(x, mu(i), sigma(i)); 
end