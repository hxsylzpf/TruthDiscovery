function diff_level_ind = func_diff_level_std(claim_std, z_std_thre)
% this function computes the diff_leve_ind according to std in claims

% num of diff levels
n_s = length(claim_std);

% if z_thre = [10 20 30]
% then the diff levels are 1: z<=10, 2: 10<z<=20, 3: 20<z<=30, 4: z>30
n_level = length(z_std_thre) + 1;

% initialize
diff_level_ind = zeros(1,n_s);

% assign level ind for each z
% the first level
idx_first = (claim_std <= z_std_thre(1));
diff_level_ind(idx_first) = 1;

% the last level
idx_last = (claim_std > z_std_thre(end));
diff_level_ind(idx_last) = n_level;

% for other levels
for i = 2:length(z_std_thre)
    cur_idx = (claim_std > z_std_thre(i-1)) & (claim_std <= z_std_thre(i));
    diff_level_ind(cur_idx) = i;
end
