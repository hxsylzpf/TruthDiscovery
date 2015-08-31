function [f, grad] = func_obj_grad_abg_global(x, n_pt, n_event, provided_label_mat, Z_est)

% x is constructed as [a b g d]
% a - x(1) to x(n_pt)
% b - x(n_pt + 1) to x(2*n_pt)
% g - x(2*n_pt + 1)
% d - x(2*n_pt + 2)

f = 0;
for j = 1:n_event
    f1 = 0;
    f2 = 0;
    for i = 1:n_pt
            f1 = f1 + provided_label_mat(i,j)*mylog(x(2*n_pt+1)*x(i)) + (1-provided_label_mat(i,j))*mylog(1- x(2*n_pt+1)*x(i));
            f2 = f2 + provided_label_mat(i,j)*mylog(x(2*n_pt+1)*x(n_pt +i)) + (1-provided_label_mat(i,j))*mylog(1 - x(2*n_pt+1)*x(n_pt+i));
    end
    f1 = f1 + mylog(x(2*n_pt + 2));
    f2 = f2 + mylog(1 - x(2*n_pt + 2));
    f = f + f1*Z_est(2,j) + f2*Z_est(1,j);
end

f = -f;

if nargout > 1
    grad_a = zeros(n_pt, 1);
    grad_b = zeros(n_pt, 1);
    grad_g = 0;
    grad_d = 0;
    
    for i = 1:n_pt
        for j = 1:n_event
            grad_a(i) = grad_a(i) + Z_est(2,j)*(provided_label_mat(i,j) - x(2*n_pt+1)*x(i))/(x(i)*(1 - x(2*n_pt+1)*x(i)));
            grad_b(i) = grad_b(i) + Z_est(1,j)*(provided_label_mat(i,j) - x(2*n_pt+1)*x(n_pt+i))/(x(n_pt+i)*(1 - x(2*n_pt+1)*x(n_pt+i)));
        end    
    end

    for j = 1:n_event
        grad_d = grad_d + Z_est(2,j)/x(2*n_pt + 2) - Z_est(1,j)/(1-x(2*n_pt + 2));
        for i = 1:n_pt
            % only has value wrt x_ij = 0; otherwise, g_ij = 1
            grad_g = grad_g + Z_est(2,j)*(provided_label_mat(i,j) - x(2*n_pt+1)*x(i))/(x(2*n_pt+1)*(1 - x(2*n_pt+1)*x(i))) ...
                                                            + Z_est(1,j)*(provided_label_mat(i,j) - x(2*n_pt+1)*x(n_pt+i))/(x(2*n_pt+1)*(1 - x(2*n_pt+1)*x(n_pt+i)));
        end
    end
    
    grad = [grad_a; grad_b; grad_g; grad_d];
    grad = -grad;
end

