function func_tbp_new_para_mapper_parfor(files, n_l, ini_pi, ini_h, ini_lambda)

%% number of files to process
n_file = numel(files);

%% for each file
% change the following to "for" or "parfor"

parfor file_idx = 1:n_file
    
    cur_file_name = files{file_idx};
    % scan the file and load data
    a = load(cur_file_name);
    X = a.X_res;
    
    [prior_mu, prior_nu] = func_prior_mu_nu(X);
    
    %% this is the mapper function
    n_pt = size(X,1);
    n_s = size(X,2);

    pi = ini_pi;
    h = ini_h;
    lambda = ini_lambda;

        %% E-step
        alpha = zeros(n_s, n_l);    
        beta = zeros(n_s, n_l);
        gamma = zeros(n_s, n_l);

        temp1 = zeros(n_s, n_l);
        temp2 = zeros(n_s, n_l);

        z_est = zeros(1, n_s);

        %% for all targets
        for j = 1:n_s       
            % for all level
            for k = 1:n_l            
                for i = 1:n_pt
                    if ~isnan(X(i,j))
                        temp1(j,k) = temp1(j,k) + lambda(i,k)*(X(i,j) - h(i,k));
                        temp2(j,k) = temp2(j,k) + lambda(i,k);
                    end
                end

                temp_mu = temp1(j,k)/temp2(j,k);
                temp_nu = temp2(j,k);

                alpha(j,k) = pi(k)*normpdf(temp_mu, prior_mu(j),sqrt(1/prior_nu(j) + 1/temp_nu));

            end        

            alpha(j,:) = alpha(j,:)/sum(alpha(j,:));
    %         show1 = alpha(j,:)

            temp3 = (prior_mu(j)*prior_nu(j) + temp1(j,:))./(prior_nu(j) + temp2(j,:));
            beta(j,:) = alpha(j,:).*temp3;
            gamma(j,:) = alpha(j,:).*(temp3.^2 + 1./(prior_nu(j) + temp2(j,:)));

            z_est(j) = sum(beta(j,:));  
        end

        %% key and val
        s1 = zeros(n_pt, n_l);
        s2 = zeros(n_pt, n_l);
        s3 = zeros(n_pt, n_l);

        % n_l vector
        s0 = sum(alpha);

        for i = 1:n_pt
            idx_i = ~isnan(X(i,:));
            for k = 1:n_l
                s1(i,k) = sum(alpha(idx_i,k));
                s2(i,k) = X(i,idx_i)*alpha(idx_i,k) - sum(beta(idx_i,k));
                s3(i,k) = X(i,idx_i).^2*alpha(idx_i,k) - 2*X(i,idx_i)*beta(idx_i,k) + sum(gamma(idx_i,k));
            end

        end

        func_tbp_mapper_save(['mapper_output_split_' num2str(n_file) '_file_' num2str(file_idx) '.mat'], n_s, s0, s1, s2, s3, z_est);

end
