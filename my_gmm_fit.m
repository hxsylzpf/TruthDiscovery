function [pi, mu, Sigma] = my_gmm_fit(data, n_l)

% this function uses kmeans to initialze pi, mu and Sigma

% data is a column vector

% n_l is the num of components
% pi is the mixture coeffient
% mu is the mixture mean
% Sigma is the mixture var

n = size(data, 1);

ini_pi = zeros(1,n_l);
ini_mu = zeros(1,n_l);
ini_Sigma = zeros(1,n_l);

if var(data) < 1e-3
   ini_pi = ones(1,n_l)/n_l;
   ini_Sigma = 1e-4*ones(1,n_l);
   
else
% k-means fit
[idx,~] = kmeans(data, n_l, 'Start', 'uniform'); % , 'Replicates',5);

for k = 1:n_l
    cur_idx = (idx == k);
    ini_pi(k) = (sum(cur_idx) + 1)/(n + n_l);
    ini_mu(k) = mean(data(cur_idx));
    ini_Sigma(k) = var(data(cur_idx));
end

end
%% GMM fit
% n_iter = 100;
% 
% gamma = zeros(n, n_l);
% 
% pi = ini_pi;
% mu = ini_mu;
% Sigma = ini_Sigma;
% 
% for kkk = 1:n_iter
%     % E-step
%     for i = 1:n
%         for k = 1:n_l
%             gamma(i,k) = pi(k)*my_normpdf(data(i),mu(k),sqrt(Sigma(k)));
%         end
%         gamma(i,:) = gamma(i,:)/sum(gamma(i,:));
%     end
%     
%     % M-step
%     new_pi = zeros(1,n_l);
%     new_mu = zeros(1,n_l);
%     new_Sigma = zeros(1,n_l);
% 
%     for k = 1:n_l
%         nk = sum(gamma(:,k));
%         
%         new_pi(k) = (nk+0.1)/(n+0.1*n_l);
%         new_mu(k) = sum(gamma(:,k).*data)/nk;
%         new_Sigma(k) = (sum(gamma(:,k).*(data - new_mu(k)).^2) + 0.1)/(nk + 0.1);
%         
%     end
%     
%     if norm(new_mu - mu) + norm(new_Sigma - Sigma) < 1e-4
%        break;
%     else
%        pi = new_pi;
%        mu = new_mu;
%        Sigma = new_Sigma;
%     end    
%     
% end

% sort by h
[~,I] = sort(abs(ini_mu));

pi = ini_pi(I);
mu = ini_mu(I);
Sigma = ini_Sigma(I);

