function record = func_solve_z_h_lambda_iter(data_table, last_val, n_l, gamma, z_est, prior_mu, prior_nu, prior_a, prior_b)

% an inner iteration for optimizing z and lambda
n_iter = 1; % 100;
% lambda_est = rand(n_pt, n_l);

n_iter_to_converge = 0;

%% find current unique pt and unique fact
cur_uni_pt = unique(data_table(:,1));
cur_uni_s = unique(data_table(:,2));

n_cur_uni_pt = length(cur_uni_pt);       
n_cur_uni_s = length(cur_uni_s);

% the max pt_id and s_id so far, they are used to allocate space
if ~isempty(last_val)
    max_n_pt = max([cur_uni_pt; last_val.max_n_pt]);
    max_n_s = max([cur_uni_s; last_val.max_n_s]);
else
    max_n_pt = max(cur_uni_pt);
    max_n_s = max(cur_uni_s);
end

%% num of current data insts
n_cur_data = size(data_table, 1);

    %% iteration
    for kkk = 1:n_iter    
    
    %% initialization - very important
    delta = zeros(max_n_pt, n_l);
    zeta = zeros(max_n_pt, n_l);
    rho = zeros(max_n_pt, n_l);
    eta = zeros(1, max_n_s);
    tau = zeros(1, max_n_s);
    
    new_z_est = zeros(size(z_est));
    new_h_est = zeros(max_n_pt, n_l);
    new_lambda_est = zeros(max_n_pt, n_l);
        
    % copy last val if not empty
    % update happens only for pts and facts in the incremental data
    if ~isempty(last_val)
        n_pt_last = last_val.max_n_pt;
        n_s_last = last_val.max_n_s;
        
        delta(1:n_pt_last,:) = last_val.delta;
        zeta(1:n_pt_last,:) = last_val.zeta;
        rho(1:n_pt_last,:) = last_val.rho;
        eta(:,1:n_s_last) = last_val.eta;
        tau(:,1:n_s_last) = last_val.tau;
        
        new_z_est(1:n_s_last) = last_val.z_est;
        new_h_est(1:n_pt_last,:) = last_val.h_est;
        new_lambda_est(1:n_pt_last,:) = last_val.lambda_est;
    end
        
    %% update delta
    % read data table rows one by one
    for nn = 1:n_cur_data
        cur_pt_id = data_table(nn,1);
        cur_s_id = data_table(nn,2);
        cur_val = data_table(nn,3);
        for k = 1:n_l
            delta(cur_pt_id,k) = delta(cur_pt_id,k) + gamma(cur_s_id,k)*(cur_val - z_est(cur_s_id));
            zeta(cur_pt_id,k) = zeta(cur_pt_id,k) + gamma(cur_s_id,k);
        end
    end  
    
    %% update h for corresponding pts; not all
    for mm = 1:n_cur_uni_pt
        cur_pt_id = cur_uni_pt(mm);
        new_h_est(cur_pt_id,:) = delta(cur_pt_id,:)./zeta(cur_pt_id,:);
    end
        
    %% update rho
    for nn = 1:n_cur_data
        cur_pt_id = data_table(nn,1);
        cur_s_id = data_table(nn,2);
        cur_val = data_table(nn,3);
        for k = 1:n_l
            rho(cur_pt_id,k) = rho(cur_pt_id,k) + 1/2*gamma(cur_s_id,k)*(cur_val - z_est(cur_s_id) - new_h_est(cur_pt_id,k))^2;
        end
    end  
    
    %% update lambda for corresponding pts; not all
    for mm = 1:n_cur_uni_pt
        cur_pt_id = cur_uni_pt(mm);
        new_lambda_est(cur_pt_id,:) = (1/2*zeta(cur_pt_id,:) + prior_a(cur_pt_id,:) - 1)./(rho(cur_pt_id,:) + prior_b(cur_pt_id,:));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %% update eta and tau
    for nn = 1:n_cur_data
        cur_pt_id = data_table(nn,1);
        cur_s_id = data_table(nn,2);
        cur_val = data_table(nn,3);
        
        for k = 1:n_l
            eta(cur_s_id) = eta(cur_s_id) + gamma(cur_s_id,k)*(new_lambda_est(cur_pt_id,k)*(cur_val - new_h_est(cur_pt_id,k)));
            tau(cur_s_id) = tau(cur_s_id) + gamma(cur_s_id,k)*new_lambda_est(cur_pt_id,k);
        end
    end
    
    %% update z for corresponding facts; not all
    for ll = 1:n_cur_uni_s
        cur_s_id = cur_uni_s(ll);
        new_z_est(cur_s_id) = (prior_nu(cur_s_id)*prior_mu(cur_s_id) + eta(cur_s_id))/(prior_nu(cur_s_id) + tau(cur_s_id));
    end
    
%     sum(new_z_est)
        
        if norm(new_z_est - z_est) < 1e-3
            n_iter_to_converge = kkk;
            break;
        else
            z_est = new_z_est;
        end
        
    end
    
%     record.z_est = round(new_z_est);
    record.z_est = new_z_est;
    record.h_est = new_h_est;
    record.lambda_est = new_lambda_est;
    
    if n_iter_to_converge ~= 0
        record.n_iter_to_converge = n_iter_to_converge;
    else
        record.n_iter_to_converge = n_iter;
    end
    
    record.delta = delta;
    record.zeta = zeta;
    record.rho = rho;
    record.eta = eta;
    record.tau = tau;
    
    