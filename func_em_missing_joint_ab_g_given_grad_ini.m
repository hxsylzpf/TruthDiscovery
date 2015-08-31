%% model parameter estimation from given data - EM
function record = func_em_missing_joint_ab_g_given_grad_ini(n_state, n_pt, n_event, provided_label_mat, loc_cover_mat, event_label, Z_est, a_est, b_est)
% provided_label_mat - n_participant * n_event
% loc_cover_mat - n_participant * n_event
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% prior s for each participant % generally not necessary - can be found from the label matrix
% s= sum(loc_cover_mat.*provided_label_mat, 2)/n_event;

    n_iter = 50;

% initialization the post prob mat
    d_est = 0.5;
    
    % solve for solutions
    lb = zeros(2*n_pt + 1, 1);
    ub = ones(2*n_pt + 1, 1);

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

    
    % grad is provided
    options = optimset('Algorithm','trust-region-reflective',...
            'Display','off','GradObj','on','UseParallel','always');
        
    % Run fmincon with starting point [¨C1,¨C1,¨C1], using the options structure:

    x0 = [a_est b_est d_est]  % Make a starting guess at the solution
    func_val = func_obj_grad_ab_g_given(x0, n_pt, n_event, provided_label_mat, loc_cover_mat, Z_est)

    tic
    [x, fval, mflag, output] = fmincon(@(x)func_obj_grad_ab_g_given(x, n_pt, n_event, provided_label_mat, loc_cover_mat, Z_est),x0,[],[],[],[],lb,ub,[],options)
    time = toc
    
    
%     % set a global w for all pt
%     x0 = [new_a_est new_b_est w0(1:n_feature)]  % Make a starting guess at the solution
%     cur_func_val = myfun_td_joint_global_w(x0, n_pt, n_event, n_feature, provided_label_mat, feature_mat, Z_est)
%     [x, fval] = fmincon(@(x)myfun_td_joint_global_w(x, n_pt, n_event, n_feature, provided_label_mat, feature_mat, Z_est), x0, [], [], [], [], lb(1:2*n_pt + n_feature), ub(1:2*n_pt + n_feature), [], options); % Call solver
    
    new_a_est = x(1:n_pt)
    new_b_est = x(n_pt+1:2*n_pt)
    new_d_est = x(2*n_pt + 1) % sum(Z_est(2,:))/n_event;
    
        if norm(new_a_est - a_est) + norm(new_b_est - b_est) + norm(new_d_est - d_est)<= 5*1e-3
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
    Z_est = zeros(2, n_event);

    for j = 1:n_event
        for i = 1:n_pt
    %         % use obs only; not 0
%              if loc_cover_mat(i,j) ~= 0
                A(j) = A(j)*(loc_cover_mat(i,j)*a_est(i))^provided_label_mat(i,j)*(1 - loc_cover_mat(i,j)*a_est(i))^(1 - provided_label_mat(i,j));
                B(j) = B(j)*(loc_cover_mat(i,j)*b_est(i))^provided_label_mat(i,j)*(1 - loc_cover_mat(i,j)*b_est(i))^(1 - provided_label_mat(i,j));
%              end
        end
        % prob of being true
        % multiply both A and B - in order to alleviate the underflow problem
         A(j) = A(j)*10^n_pt + 1e-2;
         B(j) = B(j)*10^n_pt + 1e-2;
         Z_est(2,j) = A(j)*d_est/(A(j)*d_est + B(j)*(1-d_est));
         Z_est(1,j) = 1 -Z_est(2,j);
    end
 
    disp('After E-step');
    
end

kk
Z_est

% confusion matrix
[~, confmtx] = confmat(n_state, n_event, event_label, [0 1], Z_est.')
[~, ~, ~, ~, ~, F, WAF] = binary_f_measure(confmtx)
   
    record.a_est = a_est;
    record.b_est = b_est;
    record.d_est = d_est;
    record.Z_est = Z_est;
    record.confmtx = confmtx;
    record.F = F;
    record.WAF = WAF;
    record.fval = fval;
    
