function f = func_obj_em(n_pt, n_event, provided_label_mat, loc_cover_mat, a_est, b_est, d_est, Z_est)

% x is constructed as [a b g d]
% a - x(1) to x(n_pt)
% b - x(n_pt + 1) to x(2*n_pt)
% g - x(2*n_pt + 1) to x(2*n_pt + n_event)
% d - x(2*n_pt + n_event + 1)

idx_a1 = (a_est > 1 - 1e-5);
idx_b1 = (b_est > 1 - 1e-5);
idx_a0 = (a_est < 1e-5);
idx_b0 = (b_est < 1e-5);

a_est(idx_a1) = 1 - 1e-5;
b_est(idx_b1) = 1 - 1e-5;
a_est(idx_a0) = 1e-5;
b_est(idx_b0) = 1e-5;

idx_d1 = (d_est > 1 - 1e-5);
idx_d0 = (d_est < 1e-5);
d_est(idx_d1) = 1 - 1e-5;
d_est(idx_d0) = 1e-5;

f = 0;
for j = 1:n_event
    f1 = 0;
    f2 = 0;
    for i = 1:n_pt
        if loc_cover_mat(i,j) == 1
            f1 = f1 + provided_label_mat(i,j)*log(a_est(i)) + (1-provided_label_mat(i,j))*log(1-a_est(i));
            f2 = f2 + provided_label_mat(i,j)*log(b_est(i)) + (1-provided_label_mat(i,j))*log(1-b_est(i));
        end
    end
    f1 = f1 + log(d_est);
    f2 = f2 + log(1-d_est);
    f = f + f1*Z_est(2,j) + f2*Z_est(1,j);
end

