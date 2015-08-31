function [f, grad] = myfun_obj_grad_tbpp_ori(z, s_id, n_l, data_mat, per_diff, h, lambda, w, w0, prior_mu, prior_nu)

% each z can be solved separately
% s_id is the id of this z

% n_l - num of diff levels

n_pt = size(data_mat,1);

% for all pt
temp = zeros(n_pt, n_l);

for i = 1:n_pt
    % if provided a diff level
    if ~isnan(per_diff(i,s_id))
        for k = 1:n_l
            temp(i,k) = exp(w(i,k)*z+w0(i,k));
        end
    end
end

f = 0;

%% prior on z
   f = f - 0.5*prior_nu(s_id)*(z - prior_mu(s_id))^2; % 0.5*log(prior_nu(s_id))

%% data log likelihood
       
        for i = 1:n_pt
            % if perceived diff is not nan
            if ~isnan(per_diff(i,s_id))
                % per_diff is the cur_level
                cur_level = per_diff(i,s_id);
                f = f + w(i,cur_level)*z + w0(i,cur_level) - log(sum(temp(i,:))) ...
                      - 0.5*lambda(i,cur_level)*(data_mat(i,s_id) - z - h(i,cur_level))^2;
            end
        end
    
    f = -f;
    
    f
    
    %% compute grad
if nargout > 1

    % prior
    grad = -prior_nu(s_id)*(z - prior_mu(s_id));
    
    norm_temp = prob_mat_nlz(temp, 'row');
    
    for i = 1:n_pt
        if ~isnan(per_diff(i,s_id))
            cur_level = per_diff(i,s_id);
            grad = grad - lambda(i,cur_level)*(z + h(i,cur_level) - data_mat(i,s_id)) ...
                   + w(i,cur_level) - w(i,:)*norm_temp(i,:).';
        end
    end
    
    grad = -grad;
end
