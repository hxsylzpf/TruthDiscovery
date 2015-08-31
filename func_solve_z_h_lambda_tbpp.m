function record = func_solve_z_h_lambda_tbpp(n_pt, n_s, n_l, X, zeta, z_est, prior_mu, prior_nu, prior_a, prior_b)

% an inner iteration for optimizing z and lambda
n_iter = 100;

n_iter_to_converge = 0;

for kkk = 1:n_iter    
    
    new_h_est = zeros(n_pt, n_l);
    new_lambda_est = zeros(n_pt, n_l);
    new_z_est = zeros(size(z_est));

    %% h
    for i = 1:n_pt
        idx_i = ~isnan(X(i,:));
        for k = 1:n_l
            temp1 = sum(zeta(i,idx_i,k).*(X(i,idx_i) - z_est(idx_i)));
            temp2 = sum(zeta(i,idx_i,k));
            new_h_est(i,k) = temp1/temp2;
            
            if isnan(new_h_est(i,k))
               new_h_est(i,k) = 0;
            end
            
        end
    end
    
    %% lambda
    for i = 1:n_pt
        idx_i = ~isnan(X(i,:));
        for k = 1:n_l
            temp3 = sum(zeta(i,idx_i,k));
            temp4 = sum(zeta(i,idx_i,k).*(X(i,idx_i) - z_est(idx_i) - new_h_est(i,k)).^2);
            new_lambda_est(i,k) = (0.5*temp3 + prior_a(i,k) - 1)/(0.5*temp4 + prior_b(i,k));
            
            if isnan(new_lambda_est(i,k))
               new_lambda_est(i,k) = 1e-6;
            end
            
        end
    end  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %% z
    for j = 1:n_s
        idx_j = ~isnan(X(:,j));
        temp5 = 0;
        temp6 = 0;
        for k = 1:n_l
            temp5 = temp5 + sum(zeta(idx_j,j,k).*new_lambda_est(idx_j,k).*(X(idx_j,j) - new_h_est(idx_j,k)));
            temp6 = temp6 + sum(zeta(idx_j,j,k).*new_lambda_est(idx_j,k));
        end
        new_z_est(j) = (prior_nu(j)*prior_mu(j) + temp5)/(prior_nu(j) + temp6);
    end

        if norm(new_z_est - z_est) < 1e-3
            n_iter_to_converge = kkk;
            break;
        else
            z_est = new_z_est;
        end
        
end
    
    record.z_est = new_z_est;
    record.h_est = new_h_est;
    record.lambda_est = new_lambda_est;
    
    if n_iter_to_converge ~= 0
        record.n_iter_to_converge = n_iter_to_converge;
    else
        record.n_iter_to_converge = n_iter;
    end
    