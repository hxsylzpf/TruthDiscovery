function [prior_a, prior_b] = func_tbpp_prior_lambda(n_pt, n_l, data_mat, per_diff, z)

nn = zeros(n_pt, n_l);
h = zeros(n_pt, n_l);
lambda = zeros(n_pt, n_l);

prior_a = zeros(n_pt, n_l);
prior_b = 1e-6*ones(n_pt, n_l);

for i = 1:n_pt
    for k = 1:n_l
        % nan will not be considered automatically
        idx = (per_diff(i,:) == k);
        nn(i,k) = sum(idx);
        if nn(i,k) > 0
            h(i,k) = sum(data_mat(i,idx) - z(idx))/nn(i,k);
            lambda(i,k) = (0.5*nn(i,k) + 1e-4)/(0.5*sum( (data_mat(i,idx) - z(idx) - h(i,k)).^2 ) + 1e-4);
            
            if lambda(i,k) > 100
               lambda(i,k) = 10;
            end
            
            prior_a(i,k) = nn(i,k)/2 + 1; % a_base*ones(n_pt,1)+1; % n_s_per_pt/10 + 1; % round(lambda_est12*b_base) + 1; % 10*ones(n_pt,1); %   
            prior_b(i,k) = prior_a(i,k)/lambda(i,k);
        end
    end
end

