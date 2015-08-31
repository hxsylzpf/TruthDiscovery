function record = func_avg_log_arbi_struc(n_pt, n_obj, n_claim, claim_val, claim_obj_ind, bi_pt_claim_mat)
% truth finder with arbitrary structure
% claims are mutual exclusive and each pt can only provide 1 claim for an obj

% n_obj is the number of events
% each obj can have multiple claims
% claim_obj_ind is the indicator of which object each claim is about
% e.g., claim_obj_ind = [1 1 1 2 2] means claims 1-3 are about obj 1 and 4-5 are about obj 2
% n_claim is the total claims about those objects
% in the above example, n_obj = 2 and n_claim = 5
% claim_val are the values of the claims
% e.g., if claims are about the num of people, then
% claim_val can be [10 12 14 8 9]
% bi_pt_claim_mat is in the form of n_pt*n_claim in [0,1] - indicating whether pt i provides claim j

n_iter = 100;

eps = 1e-5;

% number of claims provided by a pt % row sum
F = nansum(bi_pt_claim_mat,2);

% trust of pts
T = 0.8*ones(1,n_pt);
% belief of claims
B = 0.8*ones(1,n_claim);

last_B = inf*ones(1,n_claim);

% iterate until convergence
for kk = 1:n_iter
    
    for i = 1:n_pt
        idx_i = (bi_pt_claim_mat(i,:) == 1);
        T(i) = mylog(F(i))*sum(B(idx_i))/F(i); 
    end
    
    % normalization
    T = T/max(T);
    
    for j = 1:n_claim
        idx_j = (bi_pt_claim_mat(:,j) == 1);
        B(j) = sum(T(idx_j));
    end
    B = B/max(B);
    
    if norm(B - last_B) <= 1 - eps
        break;
    else
        last_B = B;
    end
end

% the finally estimated claim label - only one claim is true for an obj
% this label vec is [0,1]
est_claim_label = zeros(1, n_claim);
% this is the est fact value corresponds to each obj 
est_fact_per_obj = zeros(1, n_obj);

for k = 1:n_obj
    cur_claim_idx = (claim_obj_ind == k);
    cur_claim_max_conf = max(B(cur_claim_idx));
    cur_idx = (claim_obj_ind == k) & (B == cur_claim_max_conf);
    est_claim_label(cur_idx) = 1;
%     est_fact_per_obj(k) = claim_val(cur_idx);
    
    % may have multiple vals if the max_conf are the same
    est_claim_val = claim_val(cur_idx);
    n_est_claim_val = length(est_claim_val);
    % take out the first in the shuffled idx to be the final est
    idx_perm = randperm(n_est_claim_val);
    est_fact_per_obj(k) = est_claim_val(idx_perm(1));
    
end
   
    record.est_claim_label = est_claim_label;
    record.est_fact_per_obj = est_fact_per_obj;
    
    