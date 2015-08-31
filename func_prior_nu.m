function prior_nu = func_prior_nu(n_s, n_sample, samples)

% this function computes the expectation of prior_nu from samples

prior_nu = zeros(1, n_s);

for j = 1:n_s
    temp_record = zeros(n_sample, 1);
    % all the samples share the same idx_j and z_sup
    idx_j = ~isnan(samples(1,:,j));
    z_sup = sum(idx_j);
        
    for dd = 1:n_sample
        cur_x = samples(dd,:,j);
        
        % calculate median
        z_est = median(cur_x(idx_j));
        cur_var = mean((cur_x(idx_j) - z_est).^2);
        z_pre = 1/(cur_var + 0.01);
        
        % calculate 
        temp_record(dd) = z_sup.*z_pre;
    end
    prior_nu(j) = mean(temp_record);
end