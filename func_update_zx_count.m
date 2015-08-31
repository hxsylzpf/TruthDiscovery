function [n_zx_00, n_zx_01, n_zx_10, n_zx_11] = func_update_zx_count(i, n_zx_00, n_zx_01, n_zx_10, n_zx_11, last_z, cur_z, x)
% this function updates the global counts of z == 0 and z == 1

if (last_z == 0 && cur_z == 1)
   if x == 0
      n_zx_00(i) = n_zx_00(i) - 1;
      n_zx_10(i) = n_zx_10(i) + 1;
   elseif x == 1
      n_zx_01(i) = n_zx_01(i) - 1;
      n_zx_11(i) = n_zx_11(i) + 1;       
   end
   
elseif (last_z == 1 && cur_z == 0)
   if x == 0
      n_zx_00(i) = n_zx_00(i) + 1;
      n_zx_10(i) = n_zx_10(i) - 1;
   elseif x == 1
      n_zx_01(i) = n_zx_01(i) + 1;
      n_zx_11(i) = n_zx_11(i) - 1;       
   end
   
end