function record = func_tbp_new_batch(n_l, X, ini_pi, ini_h, ini_lambda, prior_mu, prior_nu, prior_a, prior_b)

%% initialization
n_iter = 100;

n_pt = size(X,1);
n_s = size(X,2);

pi = ini_pi;
h = ini_h;
lambda = ini_lambda;

for kkk = 1:n_iter
    %% E-step
    alpha = zeros(n_s, n_l);    
    beta = zeros(n_s, n_l);
    gamma = zeros(n_s, n_l);
    
    temp1 = zeros(n_s, n_l);
    temp2 = zeros(n_s, n_l);
    
    z_est = zeros(1, n_s);
    
    %% for all targets
    for j = 1:n_s
        idx_j = ~isnan(X(:,j));
        
        % for all level
        for k = 1:n_l            
%             for i = 1:n_pt
%                 if ~isnan(X(i,j))
%                     temp1(j,k) = temp1(j,k) + lambda(i,k)*(X(i,j) - h(i,k));
%                     temp2(j,k) = temp2(j,k) + lambda(i,k);
%                 end
%             end
            
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
        
        z_est(j) = sum(beta(j,:));
        
    end
    
    %% M-step
    new_pi = zeros(size(pi));
    new_h = zeros(size(h));
    new_lambda = zeros(size(lambda));
    
    for k = 1:n_l
        new_pi(k) = (sum(alpha(:,k)) + 1)/(n_s + n_l);
    end

%     for i = 1:n_pt
%         idx_i = ~isnan(X(i,:));
%         for k = 1:n_l
%             new_h(i,k) = (X(i,idx_i)*alpha(idx_i',k) - sum(beta(idx_i',k)))/sum(alpha(idx_i',k));
%             new_lambda(i,k) = (sum(alpha(idx_i',k)) + 0.5*prior_a(i,k) - 0.5)/(X(i,idx_i).^2*alpha(idx_i',k) + sum(alpha(idx_i',k))*new_h(i,k)^2 ...
%                             - 2*X(i,idx_i)*alpha(idx_i',k)*new_h(i,k) - 2*X(i,idx_i)*beta(idx_i',k) + 2*sum(beta(idx_i',k)*new_h(i,k)) + sum(gamma(idx_i',k)) ...
%                             + 0.5*prior_b(i,k)); 
%         end
%     end
    
%     temp5 = zeros(n_pt, n_l);
%     temp6 = zeros(n_pt, n_l);
%     temp7 = zeros(n_pt, n_l);
%   
%     %% for all pts - h
%     for i = 1:n_pt
%         for j = 1:n_s
%             if ~isnan(X(i,j))
%                 for k = 1:n_l
%                     temp5(i,k) = temp5(i,k) + alpha(j,k)*X(i,j) - beta(j,k);
%                     temp6(i,k) = temp6(i,k) + alpha(j,k);
%                 end
%             end
%         end
%         new_h(i,:) = temp5(i,:)./temp6(i,:);
%     end
%     
%     %% for all pts - lambda
%     for i = 1:n_pt
%         for j = 1:n_s
%             if ~isnan(X(i,j))
%                 for k = 1:n_l
%                     temp7(i,k) = temp7(i,k) + alpha(j,k)*(X(i,j) - new_h(i,k))^2 - 2*beta(j,k)*(X(i,j) - new_h(i,k)) + gamma(j,k);
%                 end
%             end
%         end
%         new_lambda(i,:) = (0.5*temp6(i,:) + prior_a(i,:) - 1)./(0.5*temp7(i,:) + prior_b(i,:));
%     end
        
s1 = zeros(n_pt, n_l);
s2 = zeros(n_pt, n_l);
s3 = zeros(n_pt, n_l);

    for i = 1:n_pt
        idx_i = ~isnan(X(i,:));
        for k = 1:n_l
            s1(i,k) = sum(alpha(idx_i,k));
            s2(i,k) = X(i,idx_i)*alpha(idx_i,k) - sum(beta(idx_i,k));
            s3(i,k) = X(i,idx_i).^2*alpha(idx_i,k) - 2*X(i,idx_i)*beta(idx_i,k) + sum(gamma(idx_i,k));
                
            new_h(i,k) = s2(i,k)/s1(i,k);
            new_lambda(i,k) = (0.5*s1(i,k) + prior_a(i,k) - 1)/(0.5*s1(i,k)*new_h(i,k)^2 - s2(i,k)*new_h(i,k) + 0.5*s3(i,k) + prior_b(i,k));
        end 
    end
    
    %% check convergence
        if norm(new_h(:) - h(:)) + norm(new_lambda(:) - lambda(:)) < 1e-3
            break;
        else
            pi = new_pi;
            h = new_h;
            lambda = new_lambda;
        end
end

record.alpha = alpha;
record.beta = beta;
record.gamma = gamma;
record.z_est = z_est;
record.pi_est = pi;
record.h_est = h;
record.lambda_est = lambda;
record.n_iter = kkk;


