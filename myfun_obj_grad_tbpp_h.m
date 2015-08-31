function [f, grad] = myfun_obj_grad_tbpp_h(h, pt_id, n_l, data_mat, per_diff, z, lambda, w, w0)

% each h can be solved separately
% h - n_l*1
% pt_id is the id of this pt

% n_l - num of diff levels

n_s = size(data_mat,2);

% for all facts
temp = zeros(n_s, n_l);

for j = 1:n_s
    % if provided a diff level
    if ~isnan(per_diff(pt_id,j))
        % note: h_k is added there
        for k = 1:n_l
            temp(j,k) = exp(w(pt_id,k)*(z(j)+h(k))+w0(pt_id,k)); 
        end
    end
end

f = 0;

%% data log likelihood
       
        for j = 1:n_s
            % if perceived diff is not nan
            if ~isnan(per_diff(pt_id,j))
                % per_diff is the cur_level
                cur_level = per_diff(pt_id,j);
                f = f + w(pt_id,cur_level)*(z(j)+h(cur_level)) + w0(pt_id,cur_level) - log(sum(temp(j,:))) ...
                      - 0.5*lambda(pt_id,cur_level)*(data_mat(pt_id,j) - z(j) - h(cur_level))^2;
            end
        end
    
    f = -f;
    
    %% compute grad
if nargout > 1

    % prior
    grad = zeros(n_l,1);
    
    norm_temp = prob_mat_nlz(temp, 'row');
    
   for j = 1:n_s
       if ~isnan(per_diff(pt_id,j))
           cur_level = per_diff(pt_id,j);
           grad(cur_level) = grad(cur_level) - lambda(pt_id,cur_level)*(z(j) + h(cur_level) - data_mat(pt_id,j)) ...
                  + w(pt_id,cur_level);
           for k = 1:n_l
              grad(k) = grad(k) - w(pt_id,k)*norm_temp(j,k);
           end
              
       end
    end
    
%     idx_s = ~isnan(per_diff(pt_id,:));
%     for k = 1:n_l
%         grad(k) = grad(k) - w(pt_id,k)*sum(norm_temp(idx_s,k));
%     end
    
    grad = -grad;
end
