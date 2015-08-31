function record = func_GLAD_mul_truth(n_state, n_pt, n_event, provided_label_mat, event_label, loc_cover_mat, Z_est)
%% model parameter estimation from given data - EM

% initialization
alpha_est = 0.5 + randn(1, n_pt);
% make alpha normalized
alpha_est = alpha_est/norm(alpha_est);
beta0_est = 1 + randn(1, n_event); % actually estimate beta0_est in the following
% % each col corresponds to an event - state 1: false; state 2: true
z_prior = 0.5*ones(n_state,n_event);
% % initialization the post prob mat
% z_post = z_prior;

n_iter = 50;

% E-step - pseudo label prob
% initialize current correct prob mat
z_post = Z_est;

for kk = 1:n_iter

%% M-step - update para est
% solve for a system of equations

x0 = [alpha_est beta0_est];  % Make a starting guess at the solution
options = optimset('Display','off', 'Algorithm', 'levenberg-marquardt', 'UseParallel', 'always'); % Option to display output % optimset('Display','iter')
[x,~] = fsolve(@(x)myfun_td(x, n_pt, n_event, provided_label_mat, loc_cover_mat, z_post), x0, options); % Call solver

new_alpha_est = x(1:n_pt);
new_beta0_est = x(n_pt + 1 : n_pt + n_event);

% make alpha normalized - otherwise, alpha and beta can be arbitrary values
new_alpha_est = new_alpha_est/norm(new_alpha_est);

    if norm(new_alpha_est - alpha_est) + norm(new_beta0_est - beta0_est) < 1e-4
        break;
    else
        alpha_est = new_alpha_est;
        beta0_est = new_beta0_est;
    end
    
%% E-step
cur_correct_prob_mat = zeros(n_pt, n_event);
z_post = zeros(2, n_event);

for j = 1:n_event
    z_post(:, j) = z_prior(:, j);
    for i = 1:n_pt
        cur_correct_prob = 1/(1+exp(-alpha_est(i)*exp(beta0_est(j))));
        cur_correct_prob_mat(i,j) = cur_correct_prob;
        % revise para only for covered events
        if loc_cover_mat(i,j) == 1
            if provided_label_mat(i,j) == 0 % first row - true label is 0 and second 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                z_post(:,j) = z_post(:,j).*[cur_correct_prob; 1 - cur_correct_prob];
    %             % revised as do nothing if label is 0
    %             disp('0');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else if provided_label_mat(i,j) == 1
                    z_post(:,j) = z_post(:,j).*[1 - cur_correct_prob; cur_correct_prob];
                end
            end
        end
    end
    % add a small positive smooting val
    z_post(:,j) = z_post(:,j)*10^n_pt + 1e-3*ones(2,1);
    % normalize corresponding col
    z_post(:,j) = z_post(:,j)/sum(z_post(:,j));
end



end

% 1/beta -> 0 - task is ambiguous; 1/beta -> inf - task is very easy
beta_est = exp(beta0_est);

alpha_est
beta_est
z_post
fval = func_obj_glad(n_pt, n_event, provided_label_mat, alpha_est, beta_est, z_post)

% confusion matrix
[~, confmtx] = confmat(2, n_event, event_label, [0 1], z_post.')
[~, ~, ~, ~, ~, F, WAF] = binary_f_measure(confmtx);

    record.confmtx = confmtx;
    record.F = F;
    record.WAF = WAF;
    record.fval = fval;
    record.est_label = (z_post > 0.5);
    