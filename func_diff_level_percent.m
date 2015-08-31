function diff_level_ind = func_diff_level_percent(z_est, z_thre)

% num of diff levels

n_s = length(z_est);

% if z_thre = [10 40]
% then the diff levels are 1: z<=10 or z>=100-10, 2: 10<z<=40 or 100-40<z<=100-10 3: 40<z<100-40
n_level = length(z_thre) + 1;

% initialize
diff_level_ind = zeros(1,n_s);

% assign level ind for each z
% the first level
idx_first = (z_est <= z_thre(1)) | (z_est >= 100 - z_thre(1));
diff_level_ind(idx_first) = 1;

% the last level
idx_last = (z_thre(end) <= z_est) & (z_est <= 100 - z_thre(end));
diff_level_ind(idx_last) = n_level;

% for other levels
for i = 2:length(z_thre)
    cur_idx = ((z_est > z_thre(i-1)) & (z_est < z_thre(i))) | ((z_est > 100 - z_thre(i)) & (z_est < 100 - z_thre(i - 1)));
    diff_level_ind(cur_idx) = i;
end