function f = myfun_td_abg(x, n_pt, n_event, provided_label_mat, Z_est)

% x is constructed as [a b g]
% a - x(1) to x(n_pt)
% b - x(n_pt + 1) to x(2*n_pt)
% g - x(2*n_pt + 1) to x(2*n_pt + n_event)
% d - x(2*n_pt + n_event + 1)

f = 0;
for j = 1:n_event
    f1 = 0;
    f2 = 0;
    for i = 1:n_pt
        f1 = f1 + provided_label_mat(i,j)*mylog(x(2*n_pt+j)*x(i)) + (1-provided_label_mat(i,j))*mylog(1- x(2*n_pt+j)*x(i));
        f2 = f2 + provided_label_mat(i,j)*mylog(x(2*n_pt+j)*x(n_pt +i)) + (1-provided_label_mat(i,j))*mylog(1 - x(2*n_pt+j)*x(n_pt+i));
    end
    f1 = f1 + mylog(x(2*n_pt + n_event + 1));
    f2 = f2 + mylog(1 - x(2*n_pt + n_event + 1));
    f = f + f1*Z_est(2,j) + f2*Z_est(1,j);
end

f = -f;
