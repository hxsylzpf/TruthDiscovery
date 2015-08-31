function result = real_to_bi_claim(pt_obj_claim_mat)
% pt_obj_claim_mat - n_pt*n_obj

n_obj = size(pt_obj_claim_mat,2);

n_claim = 0;
claim_val = [];
claim_obj_ind = [];
bi_pt_claim_mat = [];
uni_val_store = cell(n_obj, 1);

for k = 1:n_obj
    idx_k = ~isnan(pt_obj_claim_mat(:,k));
    uni_val_store{k} = unique(pt_obj_claim_mat(idx_k,k));
    cur_n_uni = length(uni_val_store{k});
    
    claim_val = [claim_val uni_val_store{k}.'];
    claim_obj_ind = [claim_obj_ind k*ones(1,cur_n_uni)];
    n_claim = n_claim + cur_n_uni;
    
    for j = 1:cur_n_uni
       cur_claim_val = uni_val_store{k}(j);
       % who claimed this val
       cur_pt_idx = (pt_obj_claim_mat(:,k) == cur_claim_val);
       bi_pt_claim_mat = [bi_pt_claim_mat cur_pt_idx];
    end    
end

result.claim_val = claim_val;
result.n_claim = n_claim;
result.n_obj = n_obj;
result.claim_obj_ind = claim_obj_ind;
result.bi_pt_claim_mat = bi_pt_claim_mat;
