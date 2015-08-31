function [n_res_per_worker, idx_retain_worker, idx_retain_res] = num_res_filter(worker_id, uni_worker_id, n_min_res)

% filter by num of responses per worker
n_uni_worker = length(uni_worker_id);
n_res_per_worker = zeros(n_uni_worker, 1);
for i = 1:n_uni_worker
    n_res_per_worker(i) = sum(strcmp(worker_id, uni_worker_id{i}));
end

% which worker should be retained
idx_retain_worker = (n_res_per_worker >= n_min_res);
retain_uni_worker_id = uni_worker_id(idx_retain_worker);

% which res should be retained
idx_retain_res = ismember(worker_id, retain_uni_worker_id);
