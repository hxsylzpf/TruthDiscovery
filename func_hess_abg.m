function hess = func_hess_abg(x, lambda, n_pt, n_event, provided_label_mat, Z_est)

% x is constructed as [a b g d]
% a - x(1) to x(n_pt)
% b - x(n_pt + 1) to x(2*n_pt)
% g - x(2*n_pt + 1) to x(2*n_pt + n_event)
% d - x(2*n_pt + n_event + 1)

hess_aa = zeros(n_pt, n_pt);
hess_bb = zeros(n_pt, n_pt);
hess_gg = zeros(n_event, n_event);
hess_dd = 0;

hess_ab = zeros(n_pt, n_pt);
hess_ag = zeros(n_pt, n_event);
hess_ad = zeros(n_pt, 1);
hess_bg = zeros(n_pt, n_event);
hess_bd = zeros(n_pt, 1);
hess_gd = zeros(n_event, 1);

for i = 1:n_pt
    for j = 1:n_event
        hess_aa(i,i) = hess_aa(i,i) + Z_est(2,j)*(2*provided_label_mat(i,j)*x(2*n_pt + j)*x(i) - provided_label_mat(i,j) - (x(2*n_pt + j)*x(i))^2)/(x(i)*(1 - x(2*n_pt + j)*x(i)))^2;
        hess_bb(i,i) = hess_bb(i,i) + Z_est(1,j)*(2*provided_label_mat(i,j)*x(2*n_pt + j)*x(n_pt+i) - provided_label_mat(i,j) - (x(2*n_pt + j)*x(n_pt+i))^2)/(x(n_pt+i)*(1 - x(2*n_pt + j)*x(n_pt+i)))^2;
        
        hess_ag(i,j) = Z_est(2,j)*(provided_label_mat(i,j) - 1)/(1 - x(2*n_pt + j)*x(i))^2;
        hess_bg(i,j) = Z_est(1,j)*(provided_label_mat(i,j) - 1)/(1 - x(2*n_pt + j)*x(n_pt+i))^2;
    end
end

for j = 1:n_event
    for i = 1:n_pt
        hess_gg(j,j) = hess_gg(j,j) + Z_est(2,j)*(2*provided_label_mat(i,j)*x(2*n_pt + j)*x(i) - provided_label_mat(i,j) - (x(2*n_pt + j)*x(i))^2)/(x(2*n_pt + j)*(1 - x(2*n_pt + j)*x(i)))^2 ...
                                                          + Z_est(1,j)*(2*provided_label_mat(i,j)*x(2*n_pt + j)*x(n_pt+i) - provided_label_mat(i,j) - (x(2*n_pt + j)*x(n_pt+i))^2)/(x(2*n_pt + j)*(1 - x(2*n_pt + j)*x(n_pt+i)))^2;
    end
end

hess = [hess_aa hess_ab hess_ag hess_ad;
               hess_ab.' hess_bb hess_bg hess_bd;
               hess_ag.' hess_bg.' hess_gg hess_gd;
               hess_ad.' hess_bd.' hess_gd.' hess_dd];
hess = - hess;

