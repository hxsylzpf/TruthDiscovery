function record = func_tbp_new_incre(n_l, X, ini_pi, ini_h, ini_lambda, prior_mu, prior_nu, prior_a, prior_b, step_exp)

%% initialization
n_pt = size(X,1);
n_s = size(X,2);

pi = ini_pi;
h = ini_h;
lambda = ini_lambda;

alpha = zeros(n_s, n_l);    
beta = zeros(n_s, n_l);
gamma = zeros(n_s, n_l);
    
temp1 = zeros(n_s, n_l);
temp2 = zeros(n_s, n_l);
   
s0 = zeros(1, n_l);
s1 = zeros(n_pt, n_l);
s2 = zeros(n_pt, n_l);
s3 = zeros(n_pt, n_l);

% num of all the claims
n_total_claim = 0;
% num of claims for each pt
n_claim = zeros(n_pt,1);

n_min = 20;

% integer idx for all the pts
idx_all_int = 1:n_pt;

% update on each new target
for j = 1:n_s
    %% E-step

    n_total_claim = n_total_claim + 1;
    
    idx_j = ~isnan(X(:,j));
    
    z_sup = sum(idx_j);
    
    % integer idx
    idx_j_int = idx_all_int(idx_j);
        
    n_claim(idx_j) = n_claim(idx_j) + 1;
     
    % for all level
    for k = 1:n_l
        temp1(j,k) = lambda(idx_j,k).'*(X(idx_j,j) - h(idx_j,k));
        temp2(j,k) = sum(lambda(idx_j,k));
            
        temp_mu = temp1(j,k)/temp2(j,k);
        temp_nu = temp2(j,k);
            
        alpha(j,k) = pi(k)*normpdf(temp_mu,prior_mu(j),sqrt(1/prior_nu(j) + 1/temp_nu));
    end    
        
        alpha(j,:) = alpha(j,:)/sum(alpha(j,:));
%         show1 = alpha(j,:)
                
        temp3 = (prior_mu(j)*prior_nu(j) + temp1(j,:))./(prior_nu(j) + temp2(j,:));
        beta(j,:) = alpha(j,:).*temp3;
        gamma(j,:) = alpha(j,:).*(temp3.^2 + 1./(prior_nu(j) + temp2(j,:)));
        
    %% update sufficient stats
    step_pi = (n_total_claim)^(-step_exp);
              
    for k = 1:n_l
        s0(k) = (1-step_pi)*s0(k) + step_pi*alpha(j,k);
        for cur_idx = 1:z_sup
            i = idx_j_int(cur_idx);
%         for i = 1:n_pt
%             if ~isnan(X(i,j))
                step = (n_claim(i))^(-step_exp);
                s1(i,k) = (1-step)*s1(i,k) + step*alpha(j,k);
                s2(i,k) = (1-step)*s2(i,k) + step*(alpha(j,k)*X(i,j) - beta(j,k));
                s3(i,k) = (1-step)*s3(i,k) + step*(alpha(j,k)*X(i,j)^2 - 2*beta(j,k)*X(i,j) + gamma(j,k));
%             end
        end
    end
            
    %% M-step
    
    if n_total_claim > n_min
        % update model para
        pi = s0;
        for cur_idx = 1:z_sup
            i = idx_j_int(cur_idx);
            if n_claim(i) > n_min
                for k = 1:n_l                
                    h(i,k) = s2(i,k)/s1(i,k);
                    lambda(i,k) = (0.5*s1(i,k) + prior_a(i,k) - 1)/(0.5*s1(i,k)*h(i,k)^2 - s2(i,k)*h(i,k) + 0.5*s3(i,k) + prior_b(i,k));
                end
            end
        end
    end
end


%% compute z
alpha = zeros(n_s, n_l);
beta = zeros(n_s, n_l);
z_est = zeros(1, n_s);

temp1 = zeros(n_s, n_l);
temp2 = zeros(n_s, n_l);

    for j = 1:n_s
        idx_j = ~isnan(X(:,j));
        % for all level
        for k = 1:n_l            
            temp1(j,k) = lambda(idx_j,k).'*(X(idx_j,j) - h(idx_j,k));
            temp2(j,k) = sum(lambda(idx_j,k));
            
            temp_mu = temp1(j,k)/temp2(j,k);
            temp_nu = temp2(j,k);
            
            alpha(j,k) = pi(k)*normpdf(temp_mu,prior_mu(j),sqrt(1/prior_nu(j) + 1/temp_nu));
        end        
        
        alpha(j,:) = alpha(j,:)/sum(alpha(j,:));
%         show1 = alpha(j,:)
                
        temp3 = (prior_mu(j)*prior_nu(j) + temp1(j,:))./(prior_nu(j) + temp2(j,:));
        beta(j,:) = alpha(j,:).*temp3;
        
        z_est(j) = sum(beta(j,:));
    end

record.alpha = alpha;
record.beta = beta;
record.gamma = gamma;
record.z_est = z_est;
record.pi_est = pi;
record.h_est = h;
record.lambda_est = lambda;

