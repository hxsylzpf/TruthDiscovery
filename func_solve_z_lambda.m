function record = func_solve_z_lambda(n_pt, n_s, n_l, X, alpha, z_est, prior_mu, prior_nu, prior_a, prior_b)

% an inner iteration for optimizing z and lambda
n_iter = 100;
new_lambda_est = zeros(n_pt, n_l);
new_z_est = zeros(size(z_est));
n_iter_to_converge = 0;

    for kkk = 1:n_iter
    
    %% lambda
    for i = 1:n_pt
        idx_i = (~isnan(X(i,:)));
        for k = 1:n_l
            temp_vec = (X(i,idx_i) - z_est(idx_i)).^2;
            new_lambda_est(i,k) = (1/2*sum(alpha(idx_i,k)) + prior_a(i,k) - 1)/...
                (1/2*temp_vec*alpha(idx_i,k) + prior_b(i,k));
        end
    end

%     sum(new_lambda_est)
    
    %% z
    for j = 1:n_s
        idx_j = (~isnan(X(:,j)));
        % dim is 1*n_l
        temp_vec_s = X(idx_j,j).'*new_lambda_est(idx_j,:);
        temp_vec_t = sum(new_lambda_est(idx_j,:), 1);
        new_z_est(j) = (prior_nu(j)*prior_mu(j) + alpha(j,:)*temp_vec_s.')/(prior_nu(j) + alpha(j,:)*temp_vec_t.');
    end
    
%     sum(new_z_est)
        
        if norm(new_z_est - z_est) < 1e-4
            n_iter_to_converge = kkk;
            break;
        else
            z_est = new_z_est;
        end
        
    end
    
%     record.z_est = round(new_z_est);
    record.z_est = new_z_est;
    record.lambda_est = new_lambda_est;
    if n_iter_to_converge ~= 0
        record.n_iter_to_converge = n_iter_to_converge;
    else
        record.n_iter_to_converge = n_iter;
    end
    