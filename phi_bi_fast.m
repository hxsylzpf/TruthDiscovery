function phi = phi_bi_fast(cur_h, cur_z, cur_x, lam_a1, lam_a0, lam_b1, lam_b0, pt_count_111, pt_count_110, pt_count_101, pt_count_100)

if cur_h == 1 && cur_z == 1 && cur_x == 1
    pt_count_111 = pt_count_111 + 1;
    pt_count_110 = pt_count_110 + 1;
    phi = (pt_count_111 + lam_a1)/(pt_count_111 + pt_count_110 + lam_a1 + lam_a0);
elseif cur_h == 1 && cur_z == 1 && cur_x == 0
    pt_count_111 = pt_count_111 + 1;
    pt_count_110 = pt_count_110 + 1;
    phi = (pt_count_110 + lam_a0)/(pt_count_111 + pt_count_110 + lam_a1 + lam_a0);
elseif cur_h == 1 && cur_z == 0 && cur_x == 1
    pt_count_101 = pt_count_101 + 1;
    pt_count_100 = pt_count_100 + 1;
    phi = (pt_count_101 + lam_b1)/(pt_count_101 + pt_count_100 + lam_b1 + lam_b0);
elseif cur_h == 1 && cur_z == 0 && cur_x == 0
    pt_count_101 = pt_count_101 + 1;
    pt_count_100 = pt_count_100 + 1;
    phi = (pt_count_100 + lam_b0)/(pt_count_101 + pt_count_100 + lam_b1 + lam_b0);
end