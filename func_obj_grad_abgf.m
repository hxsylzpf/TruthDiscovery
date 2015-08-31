function [f, grad] = func_obj_grad_abgf(x, n_pt, n_event, provided_label_mat, personal_loc_tend_mat, Z_est, eta, a_prior, b_prior, g_prior)

% x is constructed as [a b g d]
% a - x(1) to x(n_pt)
% b - x(n_pt + 1) to x(2*n_pt)
% g - x(2*n_pt + 1) to x(2*n_pt + n_event)
% d - x(2*n_pt + n_event + 1)

f = 0;
for j = 1:n_event
    f1 = 0;
    f2 = 0;
    for i = 1:n_pt
            cur_w = eta*x(2*n_pt+j) + (1-eta)*personal_loc_tend_mat(i,j);
            f1 = f1 + provided_label_mat(i,j)*log(cur_w*x(i)) + (1-provided_label_mat(i,j))*log(1- cur_w*x(i));
            f2 = f2 + provided_label_mat(i,j)*log(cur_w*x(n_pt+i)) + (1-provided_label_mat(i,j))*log(1 - cur_w*x(n_pt+i));
    end
    f1 = f1 + log(x(2*n_pt + n_event + 1));
    f2 = f2 + log(1 - x(2*n_pt + n_event + 1));
    f = f + f1*Z_est(2,j) + f2*Z_est(1,j);
    % add g_prior
    f = f + (g_prior(j,2) - 1)*log(x(2*n_pt+j)) + (g_prior(j,1) - 1)*log(1 - x(2*n_pt+j));
end

% add a_prior and b_prior
for i = 1:n_pt
    f = f + (a_prior(2) - 1)*log(x(i)) + (a_prior(1) - 1)*log(1 - x(i));
    f = f + (b_prior(2) - 1)*log(x(n_pt+i)) + (b_prior(1) - 1)*log(1 - x(n_pt+i));
end

f = -f;

if nargout > 1
    grad_a = zeros(n_pt, 1);
    grad_b = zeros(n_pt, 1);
    grad_g = zeros(n_event, 1);
    grad_d = 0;
    
    for i = 1:n_pt
        for j = 1:n_event
            cur_w = eta*x(2*n_pt+j) + (1-eta)*personal_loc_tend_mat(i,j);
            grad_a(i) = grad_a(i) + Z_est(2,j)*(provided_label_mat(i,j) - cur_w*x(i))/(x(i)*(1 - cur_w*x(i)));
            grad_b(i) = grad_b(i) + Z_est(1,j)*(provided_label_mat(i,j) - cur_w*x(n_pt+i))/(x(n_pt+i)*(1 - cur_w*x(n_pt+i)));
        end
        % add a_prior and b_prior
        grad_a(i) = grad_a(i) + (a_prior(2) - 1)/x(i) - (a_prior(1) - 1)/(1 - x(i));
        grad_b(i) = grad_b(i) + (b_prior(2) - 1)/x(n_pt+i) - (b_prior(1) - 1)/(1 - x(n_pt+i));
    end

    for j = 1:n_event
        grad_d = grad_d + Z_est(2,j)/x(2*n_pt + n_event + 1) - Z_est(1,j)/(1-x(2*n_pt + n_event + 1));
        for i = 1:n_pt
            cur_w = eta*x(2*n_pt+j) + (1-eta)*personal_loc_tend_mat(i,j);
            % only has value wrt x_ij = 0; otherwise, g_ij = 1
            grad_g(j) = grad_g(j) + Z_est(2,j)*eta*(provided_label_mat(i,j) - cur_w*x(i))/(cur_w*(1 - cur_w*x(i))) ...
                                  + Z_est(1,j)*eta*(provided_label_mat(i,j) - cur_w*x(n_pt+i))/(cur_w*(1 - cur_w*x(n_pt+i)));
        end
        % add g_prior % *x(2*n_pt+j)
        grad_g(j) = grad_g(j) + (g_prior(j,2) - 1)/x(2*n_pt+j) - (g_prior(j,1) - 1)/(1 - x(2*n_pt+j));
    end
    
    grad = [grad_a; grad_b; grad_g; grad_d];
    grad = -grad;
end
