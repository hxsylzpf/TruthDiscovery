function record = func_tbpp_no_w(n_l, data_mat, per_diff, ini_z, prior_mu, prior_nu)

% TBPP - where w does not depend on z

n_pt = size(data_mat,1);
n_s = size(data_mat,2);

h = zeros(n_pt, n_l);
lambda = zeros(n_pt, n_l);
nn = zeros(n_pt, n_l);

z = ini_z;

for iter = 1:50
% the following steps are iterated until convergence
%% solve for h and lambda given z
for i = 1:n_pt
    for k = 1:n_l
        % nan will not be considered automatically
        idx = (per_diff(i,:) == k);
        nn(i,k) = sum(idx);
        if nn(i,k) > 0
            h(i,k) = sum(data_mat(i,idx) - z(idx))/nn(i,k);
            lambda(i,k) = (0.5*nn(i,k) + 1e-4)/(0.5*sum( (data_mat(i,idx) - z(idx) - h(i,k)).^2 ) + 1e-4);
            
            if lambda(i,k) > 100
               lambda(i,k) = 10;
            end
            
        end
    end
end

%% solve for z given h and lambda
new_z = zeros(size(z));

for j = 1:n_s
    temp1 = prior_mu(j)*prior_nu(j);
    temp2 = prior_nu(j);
    
    for i = 1:n_pt
        if ~isnan(per_diff(i,j))
            cur_level = per_diff(i,j);
            temp1 = temp1 + lambda(i,cur_level)*(data_mat(i,j) - h(i,cur_level));
            temp2 = temp2 + lambda(i,cur_level);
        end
    end
    
    new_z(j) = temp1/temp2;
end

    %% check for convergence
    if norm(z - new_z) < 1e-4
       break;
    else
       z = new_z;
    end

end

record.z = z;
record.h = h;
record.lambda = lambda;
record.support = nn;
