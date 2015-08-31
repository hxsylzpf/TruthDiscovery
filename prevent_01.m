function new_val = prevent_01(val)

% this function revises 1 to 1-e, and 0 to e
% in order to prevent errors when using log

idx1 = (val>=1-1e-5);
idx0 = (val<=1e-5);

new_val = val;

new_val(idx1) = 1-1e-5;
new_val(idx0) = 1e-5;