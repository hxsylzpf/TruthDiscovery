function [f, grad] = myfun_obj_grad_tbpp_hidden_w(x, n_c, zeta, data)

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
        for k = 1:n_c
            f = f + zeta(:,j,k)*(x(k)*data(j) + x(n_c+k)) - log(sum(temp(j,:)));
        end
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
        
        for j = 1:n_data
        % sum over data
            grad_w(k) = grad_w(k) + zeta(:,j,k)*data(j);
            grad_w0(k) = grad_w0(k) + zeta(:,j,k);
        end
        
        grad_w(k) = grad_w(k) - data.'*norm_temp(:,k);
        grad_w0(k) = grad_w0(k) - sum(norm_temp(:,k));
        
        % normalize
        grad_w(k) = grad_w(k)/n_data;
        grad_w0(k) = grad_w0(k)/n_data;
        
        grad_w(k) = grad_w(k) - lam*x(k);
        grad_w0(k) = grad_w0(k) - lam*x(n_c+k);
        
    end
    
    grad = [grad_w; grad_w0];
    grad = -grad;
end
    