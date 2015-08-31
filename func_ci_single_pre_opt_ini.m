%% model parameter estimation from given data
function record = func_ci_single_pre_opt_ini(n_pt, n_s, X, ini_z, prior_mu, prior_nu, prior_a, prior_b)
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% initialization
[z_est, ~, lambda_est] = ini_z_lambda(X);
lambda_est = lambda_est.';

% z_est = ini_z;
% %% estimate lambda
% for i = 1:n_pt
%     idx_i = (X(i,:)~=0);
%     ni = sum(idx_i);
%     lambda_est(i) = (ni/2 + prior_a(i) - 1)/(1/2*sum((X(i,idx_i) - z_est(idx_i)).^2) + prior_b(i));
% end

%% optimization
    options = optimset('Algorithm','interior-point', ... % trust-region-reflective, fin-diff-grads, interior-point, lbfgs
            'Display','off','GradObj','on','UseParallel','always'); %  'Hessian', 'lbfgs', 

% Run fmincon with starting point [¨C1,¨C1,¨C1], using the options structure:
    x0 = [z_est lambda_est];  % Make a starting guess at the solution
    func_val = func_obj_grad_ci_spm(x0, n_pt, n_s, X, prior_mu, prior_nu, prior_a, prior_b);
    
    tic
    [x, fval, mflag, output] = fminunc(@(x)func_obj_grad_ci_spm(x, n_pt, n_s, X, prior_mu, prior_nu, prior_a, prior_b),x0,options)
    time = toc

    record.z_est = round(x(1:n_s));
    record.lambda_est = x(n_s+1:n_s+n_pt);
    