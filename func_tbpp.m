function record = func_tbpp(n_l, data_mat, per_diff, ini_z, prior_mu, prior_nu)

n_pt = size(data_mat,1);
n_s = size(data_mat,2);

h = zeros(n_pt, n_l);
lambda = zeros(n_pt, n_l);
w = zeros(n_pt, n_l);
w0 = zeros(n_pt, n_l);
nn = zeros(n_pt, n_l);

log_reg_acc_record = zeros(n_pt,1);

z = ini_z;

options = optimset('Algorithm','interior-point', ... % trust-region-reflective, fin-diff-grads, interior-point, lbfgs
        'Display','off','GradObj','off','Hessian', 'off', 'UseParallel','always', 'MaxIter', 500);

%% compute ini_h
for i = 1:n_pt
    for k = 1:n_l
        % nan will not be considered automatically
        idx = (per_diff(i,:) == k);
        nn(i,k) = sum(idx);
        if nn(i,k) > 0
            h(i,k) = sum(data_mat(i,idx) - z(idx))/nn(i,k);
        end
    end
end   
    
for iter = 1:50
% the following steps are iterated until convergence
%% solve for lambda given h and z
for i = 1:n_pt
    for k = 1:n_l
        % nan will not be considered automatically
        idx = (per_diff(i,:) == k);
        if nn(i,k) > 0
            lambda(i,k) = (0.5*nn(i,k) + 1e-4)/(0.5*sum( (data_mat(i,idx) - z(idx) - h(i,k)).^2 ) + 1e-4);
            
            if lambda(i,k) > 100
               lambda(i,k) = 0;
            end
            
        end
    end
end

%% solve for w given z+h
% you only need to change the input
for i = 1:n_pt
   % logistic regression
   idx = ~isnan(per_diff(i,:));
   cur_z = z(idx);
   % cur_z = abs(z(idx) - data_mat(i,idx));
   cur_per_diff = per_diff(i,idx);
   
   % copy
   cur_zh = cur_z;
   % add in bias
   for k = 1:n_l
       idx2 = (cur_per_diff == k);
       cur_zh(idx2) = cur_zh(idx2) + h(i,k);
   end
   
   %% Note: if there is no diff level k -- how about w and w0's values?
   % Note: cur_z and cur_per_diff must be col vectors
   % Note: the output cur_w and cur_w0 are also col vectors
   [cur_w, cur_w0] = func_logistic_reg_mul(n_l, cur_zh.', cur_per_diff.');
   w(i,:) = cur_w.';
   w0(i,:) = cur_w0.';
   
   target_est = func_logistic_reg_mul_pre(n_l, cur_w, cur_w0, cur_zh.');

   log_reg_acc_record(i) = sum(target_est == cur_per_diff.')/length(cur_per_diff);
   
end

%% solve for h given z, lambda and w
new_h = zeros(size(h));
% note that: each h can be solved separately, all the levels are solved together
for i = 1:n_pt
    h0 = h(i,:).';
    [h_est, fval_h, exitflag_h] = fminunc(@(h)myfun_obj_grad_tbpp_h(h, i, n_l, data_mat, per_diff, z, lambda, w, w0), h0, options);
    new_h(i,:) = h_est.';
end

%% solve for z given h, lambda and w
new_z = zeros(size(z));
      
% note that: each z can be solved separately
for j = 1:n_s
        z0 = z(j);  % Make a starting guess at the solution
        % func_val = myfun_obj_grad_em_tni(x0, i, n_fact, data_mat, sim_mat, n_claim_per_fact, pro_claim, claim_prob_mat, ind_ind, lam_r, lam_g);
        [z_est, fval_z, exitflag_z] = fminunc(@(z)myfun_obj_grad_tbpp_z(z, j, n_l, data_mat, per_diff, new_h, lambda, w, w0, prior_mu, prior_nu), z0, options);
        new_z(j) = z_est;
end

    %% check for convergence
    if norm(z - new_z) < 1e-4
       break;
    else
       h = new_h;
       z = new_z;
    end

end

record.z = z;
record.h = h;
record.lambda = lambda;
record.w = w;
record.w0 = w0;
record.log_reg_acc_record = log_reg_acc_record;
