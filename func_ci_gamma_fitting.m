function [a_est, b_est] = func_ci_gamma_fitting(X, z)

[n_pt, ~] = size(X);

a_est = zeros(n_pt, 1);
b_est = zeros(n_pt, 1);

for i = 1:n_pt
    idx_i = (~isnan(X(i,:)));
    sum(idx_i)
    ratio = X(i,idx_i)./(z(idx_i)+1e-6);
    pd = fitdist(ratio.','Gamma');
    a_est(i) = pd.a;
    b_est(i) = pd.b;
end
