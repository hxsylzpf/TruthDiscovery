function mre = my_mre(pre, gt)
% mean relative error (MRE)
% prediction, ground truth

% relative error
err = abs(pre - gt)./gt;
mre = mean(err);

