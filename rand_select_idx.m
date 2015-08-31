function out_idx = rand_select_idx(in_val, tar_val, n)

% this function randomly outputs the idx of n vals from in_val which equals tar_val
% the output is integer index

idx = find(in_val == tar_val);

% how many in_val equals tar_val
n_sat = length(idx);

idx_perm = randperm(n_sat);

if n <= n_sat
    out_idx = idx(idx_perm(1:n));
else
    out_idx = idx(idx_perm(1:n_sat));
end




