function [n_hzx_100, n_hzx_101, n_hzx_110, n_hzx_111] = func_update_hzx_count_after_z(i, n_hzx_100, n_hzx_101, n_hzx_110, n_hzx_111, ...
    last_z, cur_z, h, x)
% this function updates the global counts of hzx(i)

if (h == 1)
   if (last_z == 0) && (cur_z == 1)
       if x == 0
          n_hzx_100(i) = n_hzx_100(i) - 1;
          n_hzx_110(i) = n_hzx_110(i) + 1;
       elseif x == 1
          n_hzx_101(i) = n_hzx_101(i) - 1;
          n_hzx_111(i) = n_hzx_111(i) + 1;
       end
   elseif (last_z == 1) && (cur_z == 0)
       if x == 0
          n_hzx_100(i) = n_hzx_100(i) + 1;
          n_hzx_110(i) = n_hzx_110(i) - 1;
       elseif x == 1
          n_hzx_101(i) = n_hzx_101(i) + 1;
          n_hzx_111(i) = n_hzx_111(i) - 1;
       end
   end
   
end
