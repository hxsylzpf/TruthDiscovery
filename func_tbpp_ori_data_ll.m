function log_val = func_tbpp_ori_data_ll(n_pt, n_s, n_l, data_mat, per_diff, h, lambda, w, w0, z)

% this func computes the log data likelihood

log_val = 0;

for j = 1:n_s

%% data log likelihood
       
        for i = 1:n_pt
            % if perceived diff is not nan
            if ~isnan(per_diff(i,j))
                
                temp = zeros(1, n_l);
                for k = 1:n_l
                    temp(k) = exp(w(i,k)*z(j)+w0(i,k)); 
                end
                
                % per_diff is the cur_level
                cur_level = per_diff(i,j);
                log_val = log_val + w(i,cur_level)*z(j) + w0(i,cur_level) - log(sum(temp)) ...
                          + 0.5*log(lambda(i,cur_level)) - 0.5*lambda(i,cur_level)*(data_mat(i,j) - z(j) - h(i,cur_level))^2;
            end
        end
        
end


