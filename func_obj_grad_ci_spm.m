function [f, grad] = func_obj_grad_ci_spm(x, n_pt, n_s, X, prior_mu, prior_nu, prior_a, prior_b)

% x is constructed as [z lambda]
% z - x(1) to x(n_s)
% lambda - x(n_s + 1) to x(n_s + n_pt)

f = 0;
for j = 1:n_s
    f = f - 0.5*prior_nu(j)*(x(j) - prior_mu(j))^2;
    
    for i = 1:n_pt
        if ~isnan(X(i,j))
            f = f + 0.5*mylog(x(n_s+i)) - 0.5*x(n_s+i)*(X(i,j) - x(j))^2;
        end
    end
end

% add prior for lambda
for i = 1:n_pt
    f = f + (prior_a(i) - 1)*mylog(x(n_s+i)) - prior_b(i)*x(n_s+i);
end

f = -f;

if nargout > 1
    grad_z = zeros(n_s, 1);
    grad_lambda = zeros(n_pt, 1);
    
    for j = 1:n_s
        grad_z(j) = - prior_nu(j)*(x(j) - prior_mu(j));
        for i = 1:n_pt
            if ~isnan(X(i,j))
                grad_z(j) = grad_z(j) - x(n_s+i)*(x(j) - X(i,j));
            end
        end
    end
    
    for i = 1:n_pt
        grad_lambda(i) = (prior_a(i) - 1)/x(n_s+i) - prior_b(i);
        for j = 1:n_s
            if ~isnan(X(i,j))
                grad_lambda(i) = grad_lambda(i) + 0.5/x(n_s+i) - 0.5*(X(i,j) - x(j))^2;
            end
        end
    end
    
    grad = [grad_z; grad_lambda];
    grad = -grad;
end

