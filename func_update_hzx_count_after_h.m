function [n_hzx_100, n_hzx_101, n_hzx_110, n_hzx_111] = func_update_hzx_count_after_h(i, n_hzx_100, n_hzx_101, n_hzx_110, n_hzx_111, ...
    last_h, cur_h, z, x)
% this function updates the global counts of hzx(i)

% h:0->1
if (last_h == 0) && (cur_h == 1)
    if (z == 0) 
        if (x == 0)
            n_hzx_100(i) = n_hzx_100(i) + 1;
        elseif (x == 1)
            n_hzx_101(i) = n_hzx_101(i) + 1;
        end
    elseif (z == 1)
        if (x == 0)
            n_hzx_110(i) = n_hzx_110(i) + 1;
        elseif (x == 1)
            n_hzx_111(i) = n_hzx_111(i) + 1;
        end
    end
% h:1->0
elseif (last_h == 1) && (cur_h == 0)
    if (z == 0) 
        if (x == 0)
            n_hzx_100(i) = n_hzx_100(i) - 1;
        elseif (x == 1)
            n_hzx_101(i) = n_hzx_101(i) - 1;
        end
    elseif (z == 1)
        if (x == 0)
            n_hzx_110(i) = n_hzx_110(i) - 1;
        elseif (x == 1)
            n_hzx_111(i) = n_hzx_111(i) - 1;
        end
    end
end
   
