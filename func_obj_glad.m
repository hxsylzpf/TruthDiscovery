function f = func_obj_glad(n_pt, n_event, provided_label_mat, alpha, beta, z_est)

% Q function value
f = 0;

for i = 1:n_pt
    for j = 1:n_event
        cur_sig = sigmoid(alpha(i)*beta(j));
        if cur_sig == 1;
            cur_sig = 0.999;
        elseif cur_sig == 0
            cur_sig = 1e-5;
        end
        f = f + z_est(2,j)*(provided_label_mat(i,j)*mylog(cur_sig) + (1 - provided_label_mat(i,j))*mylog(1-cur_sig)) + ...
            z_est(1,j)*( (1-provided_label_mat(i,j))*mylog(cur_sig) + provided_label_mat(i,j)*mylog(1-cur_sig));
    end
end
        
