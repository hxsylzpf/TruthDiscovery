function record = func_truth_finder_arbi_struc_w_sim(n_pt, n_obj, n_claim, claim_val, claim_obj_ind, bi_pt_claim_mat, target)
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

gamma = 0.3;
eps = 1e-4;

% initialization
t0 = 0.8*ones(n_pt, 1);

% number of claims provided by a pt % row sum
F = nansum(bi_pt_claim_mat,2);

A = zeros(n_pt, n_claim);
B = zeros(n_claim, n_pt);

for i = 1:n_pt
    for j = 1:n_claim
       cur_claim = claim_val(j);
       if bi_pt_claim_mat(i,j) == 1
          A(i,j) = 1/F(i);
          B(j,i) = 1;
       else
          cur_obj = claim_obj_ind(j);
          % provided_info by pt i about this obj 
          cur_info_idx = (bi_pt_claim_mat(i,:) == 1) & (claim_obj_ind == cur_obj);
          cur_other_claim = claim_val(cur_info_idx);
          
          % if provides a diff claim
          if ~isempty(cur_other_claim)
               if strcmp(target, 'flight')
                  B(j,i) = exp(-abs(cur_claim - cur_other_claim)/60);
               end
          end
       end
    end
end

t = t0;

for i = 1:n_iter
    t_c = t;
    tau = - mylog(1 - t);
    sigma = B*tau;
    % confidence score of claims
    s = 1./(1 + exp(- gamma*sigma));
    % t must be smaller than 1
    t = A*s;
    
    if t.'*t_c/(norm(t)*norm(t_c)) >= 1 - eps
        break;
    end
end

% the finally estimated claim label - only one claim is true for an obj
% this label vec is [0,1]
est_claim_label = zeros(1, n_claim);
% this is the est fact value corresponds to each obj 
est_fact_per_obj = zeros(1, n_obj);

for k = 1:n_obj
    cur_claim_idx = (claim_obj_ind == k);
    cur_claim_max_conf = max(s(cur_claim_idx));
    cur_idx = (claim_obj_ind == k) & (s.' == cur_claim_max_conf);
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
    record.est_pt_rel = t;
    
    record.est_trust = tau;
    