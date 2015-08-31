function h = func_hess_ab_fixed_gf(x,lambda,n_pt,n_event,provided_label_mat,personal_loc_tend_mat,Z_est,eta,a_prior,b_prior,g_val)

    h = zeros(2*n_pt+1, 2*n_pt+1);
    
    for i = 1:n_pt
        for j = 1:n_event
            cur_w = eta*g_val(j) + (1-eta)*personal_loc_tend_mat(i,j);
            h(i) = h(i) + Z_est(2,j)*(-provided_label_mat(i,j)/x(i)^2 - (1-provided_label_mat(i,j))*cur_w^2/(1-cur_w*x(i))^2);
            h(n_pt+i) = h(n_pt+i) + Z_est(1,j)*(-provided_label_mat(i,j)/x(n_pt+i)^2 - (1-provided_label_mat(i,j))*cur_w^2/(1-cur_w*x(n_pt+i))^2);
%             grad_a(i) = grad_a(i) + Z_est(2,j)*(provided_label_mat(i,j) - cur_w*x(i))/(x(i)*(1 - cur_w*x(i)));
%             grad_b(i) = grad_b(i) + Z_est(1,j)*(provided_label_mat(i,j) - cur_w*x(n_pt+i))/(x(n_pt+i)*(1 - cur_w*x(n_pt+i)));
        end
        % add a_prior and b_prior
        h(i) = h(i) - (a_prior(2) - 1)/x(i)^2 - (a_prior(1) - 1)/(1 - x(i))^2;
        h(n_pt+i) = h(n_pt+i) - (b_prior(2) - 1)/x(n_pt+i)^2 - (b_prior(1) - 1)/(1 - x(n_pt+i))^2;
    end

    for j = 1:n_event
        h(2*n_pt+1) = h(2*n_pt+1) - Z_est(2,j)/x(2*n_pt + 1)^2 - Z_est(1,j)/(1-x(2*n_pt + 1))^2;
    end
    
    h = -h;
    