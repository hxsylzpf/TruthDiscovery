%% model parameter estimation from given data - EM with para initialization
function record = func_em_ini(n_state, n_participant, n_event, provided_label_mat, event_label, ini_Z)
% provided_label_mat - n_participant * n_event
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% prior s for each participant % generally not necessary - can be found from the label matrix
s= sum(provided_label_mat, 2)/n_event;

%% multiple initialization - and choose the one with the largest log likelihood
n_iter = 200;

    % initialization the post prob mat
    Z_est = ini_Z;
    
    a_est = 0.4 + 0.6*rand(1, n_participant);
    b_est = 0.5*rand(1, n_participant);
    d_est = 0.5;

    
for kk = 1:n_iter

%% M-step - update para est
new_a_est = ones(size(a_est));
new_b_est = ones(size(b_est));

for i = 1:n_participant
    new_a_est(i) = provided_label_mat(i,:)*Z_est(2,:).'/sum(Z_est(2,:));
    new_b_est(i) = (s(i)*n_event - provided_label_mat(i,:)*Z_est(2,:).')/(n_event - sum(Z_est(2,:)));
end
new_d_est = sum(Z_est(2,:))/n_event;

    if norm(new_a_est - a_est) + norm(new_b_est - b_est) < 1e-5
        break;
    else
        a_est = new_a_est;
        b_est = new_b_est;
        d_est = new_d_est;
    end
    
%% E-step - pseudo label prob
A = ones(n_event, 1);
B = ones(n_event, 1);

for j = 1:n_event
    for i = 1:n_participant
%         % use 1 only
%         if provided_label_mat(i,j) == 1
            A(j) = A(j)*a_est(i)^provided_label_mat(i,j)*(1 - a_est(i))^(1 - provided_label_mat(i,j));
            B(j) = B(j)*b_est(i)^provided_label_mat(i,j)*(1 - b_est(i))^(1 - provided_label_mat(i,j));
%         end
    end
    % prob of being true
    A(j) = A(j) + 1e-4; B(j) = B(j) + 1e-4;
     Z_est(2,j) = A(j)*d_est/(A(j)*d_est + B(j)*(1-d_est));
     Z_est(1,j) = 1 -Z_est(2,j);
end
 
end

% confusion matrix
[~, confmtx] = confmat(n_state, n_event, event_label, [0 1], Z_est.')
[~, ~, ~, ~, ~, F, WAF] = binary_f_measure(confmtx)

% record values
    record.a_est = a_est;
    record.b_est = b_est;
    record.d_est = d_est;
    record.Z_est = Z_est;
    record.confmtx = confmtx;
    record.F = F;
    record.WAF = WAF;
