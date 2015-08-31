function idx = comp_mat_vec(mat, vec)

% vec is a row vector - representing an inst
% mat - each row is an inst
n = size(mat, 1);

rep_vec = ones(n,1)*vec; % repmat(vec, n, 1);

% take the absolute val of diff
val = sum(abs(mat - rep_vec), 2);

idx = (val == 0);
