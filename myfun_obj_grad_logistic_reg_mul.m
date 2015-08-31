function [f, grad] = myfun_obj_grad_logistic_reg_mul(x, n_c, data, target)

lam = 1e-4;

% n_c - num of classes
% x(1:n_c) - w_k, x(n_c+1:2n_c) - w_k0

% number of data point - 1D data
n_data = size(data,1);

temp = zeros(n_data,n_c);

for k = 1:n_c
    temp(:,k) = exp(x(k)*data+x(n_c+k)); 
end

    f = 0;
    for j = 1:n_data
        % target is the index
        idx = target(j);
        f = f + x(idx)*data(j) + x(n_c+idx) - log(sum(temp(j,:)));
    end
    
    % normalize
    f = f/n_data;
    
    f = f - 0.5*lam*norm(x)^2;
    
    f = -f;
    
    %% compute grad
if nargout > 1
    
    grad_w = zeros(n_c, 1);
    grad_w0 = zeros(n_c, 1);
    
    norm_temp = prob_mat_nlz(temp, 'row');
    
    for k = 1:n_c
        
        % find all the targets which are k
        idx = (target == k);
        % sum over data
        grad_w(k) = sum(data(idx)) - data.'*norm_temp(:,k);
        grad_w0(k) = sum(idx) - sum(norm_temp(:,k));
        
        % normalize
        grad_w(k) = grad_w(k)/n_data;
        grad_w0(k) = grad_w0(k)/n_data;
        
        grad_w(k) = grad_w(k) - lam*x(k);
        grad_w0(k) = grad_w0(k) - lam*x(n_c+k);
        
    end
    
    grad = [grad_w; grad_w0];
    grad = -grad;
end
    