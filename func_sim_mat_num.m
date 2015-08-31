function [n_claim_per_fact, uni_claim, sim_mat] = func_sim_mat_num(data_mat, sigma, style)
% computing the sim_mat for numerical claims

n_s = size(data_mat, 2);

uni_claim = cell(n_s, 1);
n_claim_per_fact = zeros(n_s, 1);
sim_mat = cell(n_s, 1);

for j = 1:n_s
    idx_j = ~isnan(data_mat(:,j));
    cur_data = data_mat(idx_j,j);
    
    uni_claim{j} = unique(cur_data);
    n_claim_per_fact(j) = length(uni_claim{j});
    
    % to compute sim_mat
    if strcmp(style, 'sim_mat')
        for k = 1:n_claim_per_fact(j)
            cur_fact = uni_claim{j}(k);
            for m = 1:n_claim_per_fact(j)
                if m ~= k
                    cur_claim = uni_claim{j}(m);
                    sim_mat{j}(k,m) = exp(-abs(cur_fact - cur_claim)/sigma);
                end
            end

            % normalize
            sim_mat{j}(k,:) = sim_mat{j}(k,:)/sum(sim_mat{j}(k,:));
        end
    else
        sim_mat = [];
    end
    
end