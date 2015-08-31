function [f, grad] = func_obj_grad_glad(x, n_pt, n_event, provided_label_mat, loc_cover_mat, Z_est)

% x is constructed as [alpha beta]
% alpha - x(1) to x(n_pt)
% beta - x(n_pt + 1) to x(n_pt + n_event)

% Q function value
f = 0;

for i = 1:n_pt
    for j = 1:n_event
        if loc_cover_mat(i,j) == 1
        cur_sig = sigmoid(x(i)*x(n_pt+j));
        if cur_sig == 1;
            cur_sig = 0.999;
        elseif cur_sig == 0
            cur_sig = 1e-5;
        end
        f = f + Z_est(2,j)*(provided_label_mat(i,j)*mylog(cur_sig) + (1 - provided_label_mat(i,j))*mylog(1-cur_sig)) + ...
            Z_est(1,j)*( (1-provided_label_mat(i,j))*mylog(cur_sig) + provided_label_mat(i,j)*mylog(1-cur_sig));
        end
    end
end
 
f = -f;

if nargout > 1
    grad_a = zeros(n_pt, 1);
    grad_b = zeros(n_event, 1);
    
    for i = 1:n_pt
        for j = 1:n_event
            if loc_cover_mat(i,j) == 1
                grad_a(i) = grad_a(i) + (provided_label_mat(i,j)*Z_est(2,j) + (1 - provided_label_mat(i,j))*Z_est(1,j) - sigmoid(x(i)*x(n_pt+j)))*x(n_pt+j);
            end
        end
    end
   
    for j = 1:n_event
        for i = 1:n_pt
            if loc_cover_mat(i,j) == 1
                grad_b(j) = grad_b(j) + (provided_label_mat(i,j)*Z_est(2,j) + (1 - provided_label_mat(i,j))*Z_est(1,j) - sigmoid(x(i)*x(n_pt+j)))*x(i);
            end
        end
    end
    
    grad = [grad_a; grad_b];
    grad = -grad;
end
        

