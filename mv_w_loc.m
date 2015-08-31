function dec_val = mv_w_loc(provided_label_mat, loc_cover_mat)
% majority voting with loc_cover_mat

n_pt_report_per_evt = sum(provided_label_mat, 1);
n_pt_visit_per_evt = sum(loc_cover_mat, 1);

dec_val = n_pt_report_per_evt./n_pt_visit_per_evt;

