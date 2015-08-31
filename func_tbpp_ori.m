function record = func_tbpp_ori(n_l, data_mat, per_diff, ini_z, prior_mu, prior_nu, prior_a, prior_b)

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
        'Display','off','GradObj','on', 'Hessian', 'off', 'UseParallel','always', 'MaxIter', 200);

    
for iter = 1:50
% the following steps are iterated until convergence
%% solve for h and lambda given z
for i = 1:n_pt
    for k = 1:n_l
        % nan will not be considered automatically
        idx = (per_diff(i,:) == k);
        nn(i,k) = sum(idx);
        if nn(i,k) > 0
            h(i,k) = sum(data_mat(i,idx) - z(idx))/nn(i,k);
            lambda(i,k) = (0.5*nn(i,k) + prior_a(i,k) - 1)/(0.5*sum( (data_mat(i,idx) - z(idx) - h(i,k)).^2 ) + prior_b(i,k));
            
            if lambda(i,k) > 100
               lambda(i,k) = 10;
            elseif lambda(i,k) < 1e-8
                lambda(i,k) = 1e-8;
            elseif isnan(lambda(i,k))
                lambda(i,k) = 1e-8;
            end
            
        end
    end
end

%% solve for w given z
for i = 1:n_pt
   % logistic regression
   idx = ~isnan(per_diff(i,:));
   cur_z = z(idx);
   % cur_z = abs(z(idx) - data_mat(i,idx));
   cur_per_diff = per_diff(i,idx);
   
   %% Note: if there is no diff level k -- how about w and w0's values?
   % Note: cur_z and cur_per_diff must be col vectors
   % Note: the output cur_w and cur_w0 are also col vectors
   [cur_w, cur_w0] = func_logistic_reg_mul(n_l, cur_z.', cur_per_diff.');
   w(i,:) = cur_w.';
   w0(i,:) = cur_w0.';
   
   target_est = func_logistic_reg_mul_pre(n_l, cur_w, cur_w0, cur_z.');

   log_reg_acc_record(i) = sum(target_est == cur_per_diff.')/length(cur_per_diff);
   
end

%% solve for z given h, lambda and w
new_z = zeros(size(z));
      
% note that: each z can be solved separately
for j = 1:n_s
        idx = ~isnan(data_mat(:,j));
        if sum(idx) >= 3
            x0 = z(j);  % Make a starting guess at the solution
            [x, fval, exitflag] = fminunc(@(z)myfun_obj_grad_tbpp_ori(z, j, n_l, data_mat, per_diff, h, lambda, w, w0, prior_mu, prior_nu), x0, options);
            
            new_z(j) = x;
        else
            new_z(j) = z(j);
        end
end

    %% check for convergence
    if norm(z - new_z) < 1e-3
       break;
    else
       z = new_z;
    end
    
end

%% after convergence, compute data log likelihood
% the part related to z and the part related to lambda
data_ll = func_tbpp_ori_data_ll(n_pt, n_s, n_l, data_mat, per_diff, h, lambda, w, w0, z);

record.z = z;
record.h = h;
record.lambda = lambda;
record.w = w;
record.w0 = w0;
record.data_ll = data_ll;
record.log_reg_acc_record = log_reg_acc_record;
