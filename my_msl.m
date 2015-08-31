function msl = my_msl(pre, gt, thre)
% mean step loss (MSL)
% prediction, ground truth, thre for loss

n = length(pre);
% absolute error
err = abs(pre - gt);
% idx of err larger than thre
idx = (err > thre);
% mean step loss
msl = sum(idx)/n;

