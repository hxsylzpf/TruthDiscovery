function log_val = func_ci_multiple_bias_incmpl_logll(n_pt, n_s, n_l, X, pi, h, lambda, z)
% the value of incomplete data log likelihood

% prob
data_ll = zeros(1,n_s);

for j = 1:n_s
    for k = 1:n_l
        temp = pi(k);
        for i = 1:n_pt
            if ~isnan(X(i,j))
                % prod
                temp = temp * my_normpdf(X(i,j),z(j)+h(i,k),sqrt(1/lambda(i,k)));
            end   
        end
        data_ll(j) = data_ll(j) + temp;
    end
end

log_val = log(sum(data_ll));

