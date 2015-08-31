function func_tbp_new_para_batch_mapper_mapreduce(data, ~, intermKVStore, n_l, pi_val, h_val, lambda_val)

% note: each row is a sample
n_s = size(data,1);
n_pt = size(data,2);

% integer idx for all the pts
idx_all_int = 1:n_pt;

    %% E-step
    alpha_val = zeros(n_s, n_l);    
    beta_val = zeros(n_s, n_l);
    gamma_val = zeros(n_s, n_l);
    
    temp1 = zeros(n_s, n_l);
    temp2 = zeros(n_s, n_l);

    %% for all targets
    for j = 1:n_s
        
        % first, convert table to array
        cur_data = table2array(data(j,:));
        
        idx_j = ~isnan(cur_data);
        
        % integer idx
        idx_j_int = idx_all_int(idx_j);
        
        % compute prior        
        z_sup = sum(idx_j);
        z_est = median(cur_data(idx_j));
        cur_var = mean((cur_data(idx_j) - z_est).^2);
        z_pre = 1/(cur_var + 0.001);
        
        cur_prior_mu = z_est;
        cur_prior_nu = z_sup.*z_pre;
        
        % for all level
        for k = 1:n_l
            % only for pts who makes a claim
            % no need to loop over pts
            for cur_idx = 1:z_sup
                i = idx_j_int(cur_idx);
%             for i = 1:n_pt
%                 if ~isnan(cur_data(i))
                    temp1(j,k) = temp1(j,k) + lambda_val(i,k)*(cur_data(i) - h_val(i,k));
                    temp2(j,k) = temp2(j,k) + lambda_val(i,k);
%                 end
            end
            
            temp_mu = temp1(j,k)/temp2(j,k);
            temp_nu = temp2(j,k);
            
            alpha_val(j,k) = pi_val(k)*normpdf(temp_mu,cur_prior_mu,sqrt(1/cur_prior_nu + 1/temp_nu));

        end        
        
        alpha_val(j,:) = alpha_val(j,:)/sum(alpha_val(j,:));
                        
        temp3 = (cur_prior_mu*cur_prior_nu + temp1(j,:))./(cur_prior_nu + temp2(j,:));
        beta_val(j,:) = alpha_val(j,:).*temp3;
        gamma_val(j,:) = alpha_val(j,:).*(temp3.^2 + 1./(cur_prior_nu + temp2(j,:)));
        
    end
    
    %% key and val
    s1 = zeros(n_pt, n_l);
    s2 = zeros(n_pt, n_l);
    s3 = zeros(n_pt, n_l);
    
    % n_l vector
    s0 = sum(alpha_val);
        
    % add to an object, the key val pair
    add(intermKVStore, 's0', [s0 n_s]);
    
    for i = 1:n_pt
        cur_data = table2array(data(:,i));
        idx_i = ~isnan(cur_data);
        for k = 1:n_l
            s1(i,k) = sum(alpha_val(idx_i,k));
            s2(i,k) = alpha_val(idx_i,k).'*cur_data(idx_i) - sum(beta_val(idx_i,k));
            s3(i,k) = alpha_val(idx_i,k).'*cur_data(idx_i).^2 - 2*beta_val(idx_i,k).'*cur_data(idx_i) + sum(gamma_val(idx_i,k));    
        end
        % add to an object, the key val pair
        add(intermKVStore, num2str(i,'%04d'), [s1(i,:) s2(i,:) s3(i,:)]);
    end
