function [f, grad] = func_obj_grad_ab_g_given(x, n_pt, n_event, provided_label_mat, loc_cover_mat, Z_est)

% x is constructed as [a b g d]
% a - x(1) to x(n_pt)
% b - x(n_pt + 1) to x(2*n_pt)
% d - x(2*n_pt + 1)

f = 0;
for j = 1:n_event
    f1 = 0;
    f2 = 0;
    for i = 1:n_pt
            f1 = f1 + provided_label_mat(i,j)*mylog(loc_cover_mat(i,j)*x(i)) + (1-provided_label_mat(i,j))*mylog(1- loc_cover_mat(i,j)*x(i));
            f2 = f2 + provided_label_mat(i,j)*mylog(loc_cover_mat(i,j)*x(n_pt +i)) + (1-provided_label_mat(i,j))*mylog(1 - loc_cover_mat(i,j)*x(n_pt+i));
    end
    f1 = f1 + mylog(x(2*n_pt + 1));
    f2 = f2 + mylog(1 - x(2*n_pt + 1));
    f = f + f1*Z_est(2,j) + f2*Z_est(1,j);
end

f = -f;

if nargout > 1
    grad_a = zeros(n_pt, 1);
    grad_b = zeros(n_pt, 1);
    grad_d = 0;
    
    for i = 1:n_pt
        for j = 1:n_event
            grad_a(i) = grad_a(i) + Z_est(2,j)*(provided_label_mat(i,j) - loc_cover_mat(i,j)*x(i))/(x(i)*(1 - loc_cover_mat(i,j)*x(i)));
            grad_b(i) = grad_b(i) + Z_est(1,j)*(provided_label_mat(i,j) - loc_cover_mat(i,j)*x(n_pt+i))/(x(n_pt+i)*(1 - loc_cover_mat(i,j)*x(n_pt+i)));
        end    
    end

    for j = 1:n_event
        grad_d = grad_d + Z_est(2,j)/x(2*n_pt + 1) - Z_est(1,j)/(1-x(2*n_pt + 1));
    end
    
    grad = [grad_a; grad_b; grad_d];
    grad = -grad;
end

