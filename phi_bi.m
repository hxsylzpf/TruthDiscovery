function phi = phi_bi(pt_id, evt_id, H, Z, X, cur_h, cur_z, cur_x, lam_a1, lam_a0, lam_b1, lam_b0)

i = pt_id;
j = evt_id;

cur_HZX = [H(i,:).' Z X(i,:).'];
% for counting purpose only; first copy
cur_HZX_j = cur_HZX;
% remove the jth row - corresponding to the jth event
cur_HZX_j(j,:) = [];

if cur_h == 1 && cur_z == 1 && cur_x == 1
    pt_count_111 = sum(cur_HZX_j(:,1)==1 & cur_HZX_j(:,2)==1 & cur_HZX_j(:,3)==1); % a_i
    pt_count_110 = sum(cur_HZX_j(:,1)==1 & cur_HZX_j(:,2)==1 & cur_HZX_j(:,3)==0); % 1-a_i
    phi = (pt_count_111 + lam_a1)/(pt_count_111 + pt_count_110 + lam_a1 + lam_a0);
elseif cur_h == 1 && cur_z == 1 && cur_x == 0
    pt_count_111 = sum(cur_HZX_j(:,1)==1 & cur_HZX_j(:,2)==1 & cur_HZX_j(:,3)==1); % a_i
    pt_count_110 = sum(cur_HZX_j(:,1)==1 & cur_HZX_j(:,2)==1 & cur_HZX_j(:,3)==0); % 1-a_i
    phi = (pt_count_110 + lam_a0)/(pt_count_111 + pt_count_110 + lam_a1 + lam_a0);
elseif cur_h == 1 && cur_z == 0 && cur_x == 1
    pt_count_101 = sum(cur_HZX_j(:,1)==1 & cur_HZX_j(:,2)==0 & cur_HZX_j(:,3)==1); % b_i
    pt_count_100 = sum(cur_HZX_j(:,1)==1 & cur_HZX_j(:,2)==0 & cur_HZX_j(:,3)==0); % 1-b_i
    phi = (pt_count_101 + lam_b1)/(pt_count_101 + pt_count_100 + lam_b1 + lam_b0);
elseif cur_h == 1 && cur_z == 0 && cur_x == 0
    pt_count_101 = sum(cur_HZX_j(:,1)==1 & cur_HZX_j(:,2)==0 & cur_HZX_j(:,3)==1); % b_i
    pt_count_100 = sum(cur_HZX_j(:,1)==1 & cur_HZX_j(:,2)==0 & cur_HZX_j(:,3)==0); % 1-b_i
    phi = (pt_count_100 + lam_b0)/(pt_count_101 + pt_count_100 + lam_b1 + lam_b0);
end