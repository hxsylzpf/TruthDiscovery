%% model parameter estimation from given data - EM
function record = func_em_missing_joint_ini(n_state, n_pt, n_event, n_feature, provided_label_mat, feature_mat, event_label, loc_cover_mat, ini_Z)
% provided_label_mat - n_participant * n_event
% loc_cover_mat - n_participant * n_event
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% prior s for each participant % generally not necessary - can be found from the label matrix
% s= sum(loc_cover_mat.*provided_label_mat, 2)/n_event;

n_iter = 300;

% initialization the post prob mat
    Z_est = ini_Z;
    w0 = -0.5 + rand(1, n_pt*n_feature);
    
    a_est = 0.4 + 0.6*rand(1, n_pt);
    b_est = 0.5*rand(1, n_pt);
    d_est = 0.5;
    
    lb = [zeros(2*n_pt, 1); -inf*ones(n_pt*n_feature,1)];
    ub = [ones(2*n_pt, 1); inf*ones(n_pt*n_feature,1)];

for kk = 1:n_iter

    %% M-step - update para est
    new_a_est = ones(size(a_est));
    new_b_est = ones(size(b_est));
       
    for i = 1:n_pt
        % each pt observed sum(loc_cover_mat(i,:)) events
            new_a_est(i) = (provided_label_mat(i,:)*Z_est(2,:).')/(sum(loc_cover_mat(i,:).*Z_est(2,:)) + 1e-4);
            new_b_est(i) = (provided_label_mat(i,:)*Z_est(1,:).')/(sum(loc_cover_mat(i,:).*Z_est(1,:)) + 1e-4);
%         new_a_est(i) = (loc_cover_mat(i,:).*provided_label_mat(i,:)*Z_est(2,:).')/sum(loc_cover_mat(i,:).*Z_est(2,:));
%         new_b_est(i) = (s(i)*n_event - loc_cover_mat(i,:).*provided_label_mat(i,:)*Z_est(2,:).')/(sum(loc_cover_mat(i,:)) - sum(loc_cover_mat(i,:).*Z_est(2,:)));
    end
    new_d_est = sum(Z_est(2,:))/n_event;
    
    % solve for solutions
    options = optimset('Display','off', 'Algorithm', 'interior-point'); % Option to display output % optimset('Display','iter')
    kk
    
    % set w_i for each pt
    x0 = [new_a_est new_b_est w0]  % Make a starting guess at the solution
    cur_func_val = myfun_td_joint(x0, n_pt, n_event, n_feature, provided_label_mat, feature_mat, Z_est)
    [x, fval] = fmincon(@(x)myfun_td_joint(x, n_pt, n_event, n_feature, provided_label_mat, feature_mat, Z_est), x0, [], [], [], [], lb, ub, [], options); % Call solver

%     % set a global w for all pt
%     x0 = [new_a_est new_b_est w0(1:n_feature)]  % Make a starting guess at the solution
%     cur_func_val = myfun_td_joint_global_w(x0, n_pt, n_event, n_feature, provided_label_mat, feature_mat, Z_est)
%     [x, fval] = fmincon(@(x)myfun_td_joint_global_w(x, n_pt, n_event, n_feature, provided_label_mat, feature_mat, Z_est), x0, [], [], [], [], lb(1:2*n_pt + n_feature), ub(1:2*n_pt + n_feature), [], options); % Call solver
    
    new_a_est = x(1:n_pt);
    new_b_est = x(n_pt+1 : 2*n_pt);
    new_w_est = x(2*n_pt+1:end);
    
        if norm(new_a_est - a_est) + norm(new_b_est - b_est) + abs(new_d_est - d_est)< 1e-4
            disp('Here before break');
            break;
        else
            a_est = new_a_est;
            b_est = new_b_est;
            d_est = new_d_est;
            w0 = new_w_est;
            disp('Not break');
        end

    %% E-step - pseudo label prob
    A = ones(n_event, 1);
    B = ones(n_event, 1);

    for j = 1:n_event
        for i = 1:n_pt
    %         % use obs only; not 0
             if loc_cover_mat(i,j) ~= 0
                A(j) = A(j)*(loc_cover_mat(i,j)*a_est(i))^provided_label_mat(i,j)*(1 - loc_cover_mat(i,j)*a_est(i))^(1 - provided_label_mat(i,j));
                B(j) = B(j)*(loc_cover_mat(i,j)*b_est(i))^provided_label_mat(i,j)*(1 - loc_cover_mat(i,j)*b_est(i))^(1 - provided_label_mat(i,j));
             end
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
