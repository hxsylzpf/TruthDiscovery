function [f, grad] = func_obj_grad_ab_fixed_gf_sep(x, i, n_event, provided_label_mat, personal_loc_tend_mat, Z_est, eta, a_prior, b_prior, g_val)

% i is the id of the pt
% solve for a [a b] pair for only one pt
% x is constructed as [a b g d]
% a - x(1)
% b - x(2)

% g is fixed at prior values

f1 = 0;
f2 = 0;

for j = 1:n_event
    cur_w = eta*g_val(j) + (1-eta)*personal_loc_tend_mat(i,j);
    f1 = f1 + Z_est(2,j)*(provided_label_mat(i,j)*log(cur_w*x(1))+(1-provided_label_mat(i,j))*log(1-cur_w*x(1)));
    f2 = f2 + Z_est(1,j)*(provided_label_mat(i,j)*log(cur_w*x(2))+(1-provided_label_mat(i,j))*log(1-cur_w*x(2)));
end

    f = f1+f2;

    % add prior
    f = f + (a_prior(2) - 1)*log(x(1)) + (a_prior(1) - 1)*log(1 - x(1));
    f = f + (b_prior(2) - 1)*log(x(2)) + (b_prior(1) - 1)*log(1 - x(2));

    f = -f;

if nargout > 1
    grad_a = 0;
    grad_b = 0;

        for j = 1:n_event
            cur_w = eta*g_val(j) + (1-eta)*personal_loc_tend_mat(i,j);
            grad_a = grad_a + Z_est(2,j)*(provided_label_mat(i,j)/x(1) - cur_w*(1-provided_label_mat(i,j))/(1-cur_w*x(1)));
            grad_b = grad_b + Z_est(1,j)*(provided_label_mat(i,j)/x(2) - cur_w*(1-provided_label_mat(i,j))/(1-cur_w*x(2)));
%             grad_a(i) = grad_a(i) + Z_est(2,j)*(provided_label_mat(i,j) - cur_w*x(i))/(x(i)*(1 - cur_w*x(i)));
%             grad_b(i) = grad_b(i) + Z_est(1,j)*(provided_label_mat(i,j) - cur_w*x(n_pt+i))/(x(n_pt+i)*(1 - cur_w*x(n_pt+i)));
        end
        % add a_prior and b_prior
        grad_a = grad_a + (a_prior(2) - 1)/x(1) - (a_prior(1) - 1)/(1 - x(1));
        grad_b = grad_b + (b_prior(2) - 1)/x(2) - (b_prior(1) - 1)/(1 - x(2));

    
    grad = [grad_a; grad_b];
    grad = -grad;
end
