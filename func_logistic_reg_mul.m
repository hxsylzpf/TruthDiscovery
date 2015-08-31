function [w, w0] = func_logistic_reg_mul(n_c, data, target)

% multiclass logistic regression
% n_c is the number of classes
% target is the label

% relation: p(target = k|data) = exp(wk*data+wk0)/sum_k' exp(wk'*data+wk'0)

    options = optimset('Algorithm','interior-point', ... % trust-region-reflective, fin-diff-grads, interior-point, lbfgs
           'Display','off','GradObj','on', 'Hessian', 'off', 'UseParallel','always', 'MaxIter', 300);
    
        x0 = rand(2*n_c,1);  % Make a starting guess at the solution
        % func_val = myfun_obj_grad_em_tni(x0, i, n_fact, data_mat, sim_mat, n_claim_per_fact, pro_claim, claim_prob_mat, ind_ind, lam_r, lam_g);
        [x, fval, exitflag] = fminunc(@(x)myfun_obj_grad_logistic_reg_mul(x, n_c, data, target), x0, options);
        w = x(1:n_c);
        w0 = x(n_c+1:2*n_c);

%         fval
%         exitflag