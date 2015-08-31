%% model parameter estimation from given data - EM
function record = func_em_missing_joint_ab_fixed_gf_grad_sep_ini(n_state, n_pt, n_event, provided_label_mat, personal_loc_tend_mat, ...
         event_label, Z_est, a_est, b_est, eta, a_prior, b_prior, g_val)
% solve para for each pt separately!!!     
% provided_label_mat - n_pt * n_evt
% personal_loc_tend_mat - n_pt * n_evt
% loc_cover_mat - n_pt * n_evt
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% prior s for each participant % generally not necessary - can be found from the label matrix
    n_iter = 50;
    
% initializing the post prob mat
    d_est = 0.5;
    
for kk = 1:n_iter

    %% M-step - update para est
   
%     % solve for solutions
%     options = optimset('Display','off','Algorithm','interior-point'); % Option to display output % optimset('Display','iter')
%     kk
%     
%     % set w_i for each pt
%     x0 = [a_est b_est g_est d_est]  % Make a starting guess at the solution
%     func_val = myfun_td_abg(x0, n_pt, n_event, provided_label_mat, Z_est)
%     [x, fval, exitflag] = fmincon(@(x)myfun_td_abg(x, n_pt, n_event, provided_label_mat, Z_est), x0, [], [], [], [], lb, ub, [], options); % Call solver

%%    % grad is provided
%     options = optimset('Algorithm','interior-point', ... % trust-region-reflective, fin-diff-grads, interior-point, lbfgs
%             'Display','off','GradObj','on', 'MaxIter', 300, 'TolFun', 1e-3, 'TolX', 1e-5, 'UseParallel','always'); % 'Hessian', 'lbfgs', 'UseParallel','always',
    
%     options = optimoptions(@fmincon,'Algorithm','interior-point', ... % trust-region-reflective, fin-diff-grads, interior-point, lbfgs
%             'Display','off','GradObj','on', 'MaxIter', 300, 'TolFun', 1e-3, 'TolX', 1e-5, 'UseParallel','always', ...
%             'Hessian','user-supplied',...
%             'HessFcn',f_hess);
%     
%     % Run fmincon with starting point [¨C1,¨C1,¨C1], using the options structure:
% 
%     x0 = [a_est b_est d_est]  % Make a starting guess at the solution
%     func_val = func_obj_grad_ab_fixed_gf(x0, n_pt, n_event, provided_label_mat, personal_loc_tend_mat, Z_est, eta, a_prior, b_prior, g_val)
% 
%     % hess = func_hess_ab_fixed_gf(x,n_pt,n_event,provided_label_mat,personal_loc_tend_mat,Z_est,eta,a_prior,b_prior,g_val)
%     
%     tic
%     [x, fval, mflag, output] = fmincon(@(x)func_obj_grad_ab_fixed_gf(x, n_pt, n_event, provided_label_mat, personal_loc_tend_mat, ...
%         Z_est, eta, a_prior, b_prior, g_val),x0,[],[],[],[],lb,ub,[],options)
%     time = toc
%     
%     
% %     % set a global w for all pt
% %     x0 = [new_a_est new_b_est w0(1:n_feature)]  % Make a starting guess at the solution
% %     cur_func_val = myfun_td_joint_global_w(x0, n_pt, n_event, n_feature, provided_label_mat, feature_mat, Z_est)
% %     [x, fval] = fmincon(@(x)myfun_td_joint_global_w(x, n_pt, n_event, n_feature, provided_label_mat, feature_mat, Z_est), x0, [], [], [], [], lb(1:2*n_pt + n_feature), ub(1:2*n_pt + n_feature), [], options); % Call solver
%     
%     new_a_est = x(1:n_pt)
%     new_b_est = x(n_pt+1:2*n_pt)
%     new_d_est = x(2*n_pt+1) % sum(Z_est(2,:))/n_event;

%% separate solvers
options = optimset('Algorithm','interior-point', ... % trust-region-reflective, fin-diff-grads, interior-point, lbfgs
         'Display','off','GradObj','on', 'MaxIter', 300, 'TolFun', 1e-5, 'TolX', 1e-6, 'UseParallel','always'); % 'Hessian', 'lbfgs', 'UseParallel','always',

% solve for solutions
lb = zeros(2,1);
ub = ones(2,1);

new_a_est = zeros(size(a_est));
new_b_est = zeros(size(b_est));

for i = 1:n_pt
   x0 = [a_est(i) b_est(i)];  % Make a starting guess at the solution
   func_val = func_obj_grad_ab_fixed_gf_sep(x0, i, n_event, provided_label_mat, personal_loc_tend_mat, Z_est, eta, a_prior, b_prior, g_val);
   [x, fval, mflag, output] = fmincon(@(x)func_obj_grad_ab_fixed_gf_sep(x, i, n_event, provided_label_mat, personal_loc_tend_mat, ...
        Z_est, eta, a_prior, b_prior, g_val),x0,[],[],[],[],lb,ub,[],options);
   new_a_est(i) = x(1);
   new_b_est(i) = x(2);
end

new_d_est = sum(Z_est(2,:))/n_event;

        if norm(new_a_est - a_est) + norm(new_b_est - b_est)<= 1e-3
            disp('Here before break');
            break;
        else
            a_est = new_a_est;
            b_est = new_b_est;
            d_est = new_d_est;
            disp('Not break');
        end

    %% E-step - pseudo label prob
    A = ones(n_event, 1);
    B = ones(n_event, 1);

    for j = 1:n_event
        for i = 1:n_pt
            cur_w = eta*g_val(j) + (1-eta)*personal_loc_tend_mat(i,j);
    %       % use obs only; not 0
%           % if loc_cover_mat(i,j) ~= 0
            A(j) = A(j)*(cur_w*a_est(i))^provided_label_mat(i,j)*(1 - cur_w*a_est(i))^(1 - provided_label_mat(i,j));
            B(j) = B(j)*(cur_w*b_est(i))^provided_label_mat(i,j)*(1 - cur_w*b_est(i))^(1 - provided_label_mat(i,j));
%           % end
            % multiply both A and B - in order to alleviate the underflow problem
            A(j) = A(j)*10;
            B(j) = B(j)*10;
        end
        % prob of being true
         A(j) = A(j) + 1e-2;
         B(j) = B(j) + 1e-2;
         Z_est(2,j) = A(j)*d_est/(A(j)*d_est + B(j)*(1-d_est));
         Z_est(1,j) = 1 - Z_est(2,j);
    end
 
    disp('After E-step');
    
end

kk
Z_est

% confusion matrix
[~, confmtx] = confmat(n_state, n_event, event_label, [0 1], Z_est.')
[~, ~, ~, pre1, rec1, F, WAF] = binary_f_measure(confmtx)
[tnr, tpr, acc] = binary_acc(confmtx);
   
    record.a_est = a_est;
    record.b_est = b_est;
    record.d_est = d_est;
    record.Z_est = Z_est;
    record.confmtx = confmtx;
    record.F = F;
    record.WAF = WAF;
    record.fval = fval;
    record.pre = pre1;
    record.rec = rec1;
    record.tnr = tnr;
    record.tpr = tpr;
    record.acc = acc;
    
