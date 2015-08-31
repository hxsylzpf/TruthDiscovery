function [n_h0, n_h1] = func_update_h_count(j, n_h0, n_h1, last_h, cur_h)
% this function updates the global counts of h(j) == 0 and h(j) == 1

if (last_h == 0 && cur_h == 1)
   % overall count
   n_h0(j) = n_h0(j) - 1;
   n_h1(j) = n_h1(j) + 1;
   
elseif (last_h == 1 && cur_h == 0)
   % overall count
   n_h0(j) = n_h0(j) + 1;
   n_h1(j) = n_h1(j) - 1;
   
end
