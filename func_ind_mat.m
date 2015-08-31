function ind_mat = func_ind_mat(n_pt, n_s, style, pt_claim_rate, n_pt_per_target)

   ind_mat = zeros(n_pt, n_s);
   
   if strcmp(style, 'pt_claim_rate')
       
       if mean(pt_claim_rate) == 1
          ind_mat = ones(n_pt, n_s);
       else
           for i = 1:n_pt
               cur_rate = pt_claim_rate(i);
               rng(i);
               ind_mat(i,:) = discreternd([1-cur_rate cur_rate],n_s) - 1; 
           end
       end
   
   elseif strcmp(style, 'n_pt_per_target')
       for j = 1:n_s
           rng(j);
           idx_perm = randperm(n_pt);
           idx_sel = idx_perm(1:n_pt_per_target);
           ind_mat(idx_sel, j) = 1;
       end
   end
   