function [ent_per_worker, idx_retain_worker, idx_retain_res, record] = entropy_filter(worker_id, uni_worker_id, X_res, min_ent)

% filter by entropy per worker
n_uni_worker = length(uni_worker_id);
ent_per_worker = zeros(n_uni_worker, 1);
record = struct;

for i = 1:n_uni_worker
   idx_i = ~isnan(X_res(i,:));
   
   % response
   cur_res = X_res(i,idx_i);
   cur_uni_res = unique(cur_res);
   record(i).uni_res = cur_uni_res;
   record(i).num_res = mycount_unique(cur_res, cur_uni_res);
   % entropy of response
   record(i).prob_res = record(i).num_res/sum(record(i).num_res);
   record(i).ent_res = myentropy_pmf(record(i).prob_res);
   ent_per_worker(i) = record(i).ent_res;
end

% which worker should be retained
idx_retain_worker = (ent_per_worker >= min_ent);
retain_uni_worker_id = uni_worker_id(idx_retain_worker);

% which res should be retained
idx_retain_res = ismember(worker_id, retain_uni_worker_id);
