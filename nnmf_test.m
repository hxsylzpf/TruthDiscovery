% NNMF test
k = 3;
option.distance = 'kl';
% X = A*Y
n_entry = sum(sum(provided_label_mat))

idx = (provided_label_mat == 0);
provided_label_mat_w_missing = provided_label_mat;
provided_label_mat_w_missing(idx) = NaN

[A,Y,numIter,tElapsed,finalResidual] = wnmfrule(provided_label_mat,k, option)

n_entry = sum(sum(provided_label_mat))

est = A*Y;
err = sum(est(logical(provided_label_mat))) - n_entry

% est for missing entries

est_missing = est(~logical(provided_label_mat))
mean_est_missing = mean(est_missing)

% estimated loc_cover_mat
re_loc_cover_mat = ones(n_pt, n_event);
idx1 = (est < mean_est_missing);
re_loc_cover_mat(idx1) = 0
% if already in provided_label_mat
re_loc_cover_mat(~idx) = 1

re_err = sum(sum(abs(loc_cover_mat - re_loc_cover_mat)))

subplot(2,1,1)
spy(loc_cover_mat)
title('Loc cover mat')
subplot(2,1,2)
spy(re_loc_cover_mat)
title('Reloc cover mat')

