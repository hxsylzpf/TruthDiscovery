%% model parameter estimation from given data
function record = func_ci_ratio_opt_ini(n_pt, n_s, X, ini_z, ini_a, ini_b)
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% optimization
    options = optimset('Algorithm','interior-point', ... % trust-region-reflective, fin-diff-grads, interior-point, lbfgs
            'Display','off','UseParallel','always'); % 'GradObj','on', % 'Hessian', 'lbfgs', 

% Run fmincon with starting point [¨C1,¨C1,¨C1], using the options structure:
    x0 = [ini_z ini_a ini_b];  % Make a starting guess at the solution
    func_val = func_obj_ci_ratio(x0, n_pt, n_s, X);
    lb = zeros(1, n_s + 2*n_pt);
    
    tic
    [x, fval, mflag, output] = fmincon(@(x)func_obj_ci_ratio(x, n_pt, n_s, X),x0,[],[],[],[],lb,[],[],options)
    time = toc

    record.z_est = x(1:n_s);
    record.a_est = x(n_s+1:n_s+n_pt);
    record.b_est = x(n_s+n_pt+1:n_s+2*n_pt);
    