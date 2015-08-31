function record = func_incre_tbp_new(n_l, X, ini_pi, ini_h, ini_lambda, prior_mu, prior_nu, prior_a, prior_b)

%% initialization
n_pt = size(X,1);
n_s = size(X,2);

pi = ini_pi;
h = ini_h;
lambda = ini_lambda;

n_claim = 0;
n_claim_pt = zeros(n_pt,1);
n_claim_s = zeros(n_s,1);

s1 = zeros(1, n_l);
s2 = zeros(n_pt, n_l);
s3 = zeros(n_pt, n_l);
s4 = zeros(n_pt, n_l);
s5 = zeros(n_pt, n_l);
s6 = zeros(n_pt, n_l);
s7 = zeros(n_pt, n_l);
s8 = zeros(1, n_s);

% scan through all the data once
for j = 1:n_s
    for i = 1:n_pt
        if ~isnan(X(i,j))
            %% update count
            n_claim = n_claim + 1;
            n_claim_pt(i) = n_claim_pt(i) + 1;
            n_claim_s(j) = n_claim_s(j) + 1;
            
            step = (n_claim)^(-0.6);
            step_pt = (n_claim_pt(i))^(-0.6);
%             step_s = (n_claim_s(j))^(-0.6);
            
            %% E-step
            %% compute sufficient stats
            temp1 = zeros(1,n_l);
            cur_mu = zeros(1,n_l);
            cur_nu = zeros(1,n_l);
            for k = 1:n_l
                temp1(k) = pi(k)*normpdf(X(i,j),prior_mu(j)+h(i,k),sqrt(1/prior_nu(j)+1/lambda(i,k)));
                cur_mu(k) = (prior_mu(j)*prior_nu(j)+lambda(i,k)*(X(i,j)-h(i,k)))/(prior_nu(j)+lambda(i,k));
                cur_nu(k) = prior_nu(j)+lambda(i,k);
            end
            
            cur_alpha = temp1/sum(temp1);
            cur_beta = cur_alpha.*cur_mu;
            cur_gamma = cur_alpha.*(cur_mu.^2 + 1./cur_nu);
            
            %% update model paras
            for k = 1:n_l
                s1(k) = (1-step)*s1(k) + step*cur_alpha(k);
                s2(i,k) = (1-step_pt)*s2(i,k) + step_pt*cur_alpha(k);
                s3(i,k) = (1-step_pt)*s3(i,k) + step_pt*cur_beta(k);
                s4(i,k) = (1-step_pt)*s4(i,k) + step_pt*cur_gamma(k);
                s5(i,k) = (1-step_pt)*s5(i,k) + step_pt*cur_alpha(k)*X(i,j);
                s6(i,k) = (1-step_pt)*s6(i,k) + step_pt*cur_alpha(k)*X(i,j)^2;
                s7(i,k) = (1-step_pt)*s7(i,k) + step_pt*cur_beta(k)*X(i,j);
            end
            
%             s8(j) = (1-step_s)*s8(j) + step_s*sum(cur_beta);
    
            %% M-step - only perform when n_claim_pt > 10
            if n_claim_pt(i) > 10
               pi = s1;
               
               for k = 1:n_l
                   h(i,k) = (s5(i,k) - s3(i,k))/s2(i,k);
                   lambda(i,k) = (0.5*s2(i,k) + prior_a(i,k) - 1)/(0.5*s2(i,k)*h(i,k)^2 + 0.5*s4(i,k) + 0.5*s6(i,k) ...
                               + s3(i,k)*h(i,k) - s5(i,k)*h(i,k) - s7(i,k) + prior_b(i,k));
               end
                
            end
            
        end
    end
end

%% estimate z after updating model paras

    log_alpha = zeros(n_s, n_l);
    alpha = zeros(n_s, n_l);    
    beta = zeros(n_s, n_l);
    
    temp1 = zeros(n_s, n_l);
    temp2 = zeros(n_s, n_l);
    
    z = zeros(1, n_s);
    
    %% for all targets
    for j = 1:n_s       
        % for all level
        for k = 1:n_l
            log_alpha(j,k) = log(pi(k));
            temp1(j,k) = prior_mu(j)*prior_nu(j);
            temp2(j,k) = prior_nu(j);
            
            for i = 1:n_pt
                if ~isnan(X(i,j))
                   log_alpha(j,k) = log_alpha(j,k) + log(normpdf(X(i,j),prior_mu(j)+h(i,k),sqrt(1/prior_nu(j) + 1/lambda(i,k))));  
                   temp1(j,k) = temp1(j,k) + lambda(i,k)*(X(i,j) - h(i,k));
                   temp2(j,k) = temp2(j,k) + lambda(i,k);                   
                end
            end
            
        end
        
    end
    
    %% for all targets
    for j = 1:n_s
        temp3 = exp(log_alpha(j,:));
        alpha(j,:) = temp3/sum(temp3);
        
        temp4 = temp1(j,:)./temp2(j,:);
        beta(j,:) = alpha(j,:).*temp4;
        
        z(j) = sum(beta(j,:));
    end

record.z = z;
record.pi = pi;
record.h = h;
record.lambda = lambda;
record.n_claim = n_claim;
record.n_claim_pt = n_claim_pt;
record.n_claim_s = n_claim_s;
