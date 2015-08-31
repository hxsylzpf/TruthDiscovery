%% model parameter estimation from given data - EM
function record = func_em(n_state, n_participant, n_event, provided_label_mat, event_label)
% provided_label_mat - n_participant * n_event
% parameter initialization is very important for EM - bad initialization
% can result in bad estimation

%% prior s for each participant % generally not necessary - can be found from the label matrix
s= sum(provided_label_mat, 2)/n_event;

%% multiple initialization - and choose the one with the largest log likelihood
n_ini = 50;
n_iter = 200;
% record log likelihood
LL = zeros(n_ini, 1);
LL_best = -inf;

for pp = 1:n_ini
    % initialization
    a_est = 0.2 + 0.7*rand(1, n_participant);
    b_est = 0.8*rand(1, n_participant);
    d_est = 0.5;
    % each col corresponds to an event - state 1: false; state 2: true
    Z_prior = 0.5*ones(n_state,n_event);
    % initialization the post prob mat
    Z_est = Z_prior;

for kk = 1:n_iter
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

end

% confusion matrix
[~, confmtx] = confmat(n_state, n_event, event_label, [0 1], Z_est.')
[~, ~, ~, ~, ~, F, WAF] = binary_f_measure(confmtx)


% current data log likelihood
for i = 1:n_participant
    for j = 1:n_event
        if provided_label_mat(i,j) == 1
            LL(pp) = LL(pp) + log(d_est*a_est(i) + (1 - d_est)*b_est(i));
        elseif provided_label_mat(i,j) == 0
            LL(pp) = LL(pp) + log(d_est*(1-a_est(i)) + (1 - d_est)* (1 - b_est(i)));
        end
    end
end

% record the best config
if LL(pp) >= LL_best
    % update the best LL
    LL_best = LL(pp);
    
    record.a_est = a_est;
    record.b_est = b_est;
    record.d_est = d_est;
    record.Z_est = Z_est;
    record.confmtx = confmtx;
    record.F = F;
    record.WAF = WAF;
end

% % record current estimations
% record{pp,1} = a_est;
% record{pp,2} = b_est;
% record{pp,3} = d_est;
% record{pp,4} = Z_est;
% record{pp,5} = confmtx;
% record{pp,6} = [F wei_ave_F];

end