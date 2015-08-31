function new_data_mat = func_sub_sampling(data_mat, n_pt_per_task)

% this function sub sample a given data_mat with n_pt_per_task for each task
% rnd_seed is the seed used for the rand function

n_s = size(data_mat,2);

new_data_mat = nan(size(data_mat));

%% sub sampling
for j = 1:n_s
   cur_X = data_mat(:,j);
   
   % find not nan
   idx = find(isnan(cur_X) == 0);
   
   n_idx = length(idx);
   
   % sub sample if more than desired pts exist
   if n_idx > n_pt_per_task
      idx_perm = idx(randperm(n_idx));
      idx_sel = idx_perm(1:n_pt_per_task);
      % copy sel
      new_data_mat(idx_sel,j) = data_mat(idx_sel,j);
   else
      new_data_mat(:,j) = data_mat(:,j);
   end
end


