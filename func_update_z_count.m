function [n_z0, n_z1] = func_update_z_count(n_z0, n_z1, last_z, cur_z)
% this function updates the global counts of z == 0 and z == 1

if (last_z == 0 && cur_z == 1)
   % overall count
   n_z0 = n_z0 - 1;
   n_z1 = n_z1 + 1;
elseif (last_z == 1 && cur_z == 0)
   % overall count
   n_z0 = n_z0 + 1;
   n_z1 = n_z1 - 1;
end