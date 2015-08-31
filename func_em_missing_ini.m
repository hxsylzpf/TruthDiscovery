%% model parameter estimation from given data - EM
function record = func_em_missing_ini(n_state, n_pt, n_event, provided_label_mat, event_label, loc_cover_mat, ini_Z)
% provided_label_mat - n_participant * n_event
% loc_cover_mat - n_participant * n_event
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% prior s for each participant % generally not necessary - can be found from the label matrix
% s= sum(loc_cover_mat.*provided_label_mat, 2)/n_event;

n_iter = 100;

% initialization the post prob mat
    Z_est = ini_Z;

    a_est = 0.4 + 0.6*rand(1, n_pt);
    b_est = 0.5*rand(1, n_pt);
    d_est = 0.5;

for kk = 1:n_iter

    %% M-step - update para est
    new_a_est = ones(size(a_est));
    new_b_est = ones(size(b_est));

    for i = 1:n_pt
        % each pt observed sum(loc_cover_mat(i,:)) events
            new_a_est(i) = (provided_label_mat(i,:)*Z_est(2,:).')/(loc_cover_mat(i,:)*Z_est(2,:).'+1e-2);
            new_b_est(i) = (provided_label_mat(i,:)*Z_est(1,:).')/(loc_cover_mat(i,:)*Z_est(1,:).'+1e-2);
%         new_a_est(i) = (loc_cover_mat(i,:).*provided_label_mat(i,:)*Z_est(2,:).')/sum(loc_cover_mat(i,:).*Z_est(2,:));
%         new_b_est(i) = (s(i)*n_event - loc_cover_mat(i,:).*provided_label_mat(i,:)*Z_est(2,:).')/(sum(loc_cover_mat(i,:)) - sum(loc_cover_mat(i,:).*Z_est(2,:)));
    end
    new_d_est = sum(Z_est(2,:))/n_event;

        if norm(new_a_est - a_est) + norm(new_b_est - b_est) + abs(new_d_est - d_est)< 1e-4
            break;
        else
            a_est = prevent_01(new_a_est);
            b_est = prevent_01(new_b_est);
            d_est = prevent_01(new_d_est);
        end

    %% E-step - pseudo label prob
    A = ones(n_event, 1);
    B = ones(n_event, 1);
    Z_est = zeros(2, n_event);

    for j = 1:n_event
        for i = 1:n_pt
    %         % use obs only; not 0
             if loc_cover_mat(i,j) == 1
                A(j) = A(j)*(a_est(i))^provided_label_mat(i,j)*(1 - a_est(i))^(1 - provided_label_mat(i,j));
                B(j) = B(j)*(b_est(i))^provided_label_mat(i,j)*(1 - b_est(i))^(1 - provided_label_mat(i,j));
             end
        end
        % prob of being true
        % multiply both A and B - in order to alleviate the underflow problem
         A(j) = A(j)+1e-2;
         B(j) = B(j)+1e-2;
         Z_est(2,j) = A(j)*d_est/(A(j)*d_est + B(j)*(1-d_est));
         Z_est(1,j) = 1-Z_est(2,j);
    end
 
end

kk

idx_valid = ~isnan(Z_est(1,:));

if sum(idx_valid == 0)
   Z_est = rand(2, n_event);
   Z_est = prob_mat_nlz(Z_est, 'col');
end


% value of Q function
fval = func_obj_em(n_pt, n_event, provided_label_mat, loc_cover_mat, a_est, b_est, d_est, Z_est)

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
    