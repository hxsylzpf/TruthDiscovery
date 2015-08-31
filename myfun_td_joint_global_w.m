function f = myfun_td_joint_global_w(x, n_pt, n_event, n_feature, provided_label_mat, feature_mat, Z_est)

% x is constructed as [alpha beta w_1 w_2 ... w_n_pt]
% alpha [length - n_pt] - x(1) to x(n_pt)
% beta [length - n_pt] - x(n_pt + 1) to x(2*n_pt)
% w [length - n_feature] - x(2*n_pt + 1) to x(2*n_pt + n_feature)

% construct the loc_cover_mat - first copy ones
loc_cover_mat = provided_label_mat;

% w is global
cur_w = x(2*n_pt + 1 : 2*n_pt + n_feature);
% then construct prob est according to logistic regression
for ii = 1:n_pt
    for jj = 1:n_event
        if loc_cover_mat(ii, jj) == 0 % not a revealed loc
            loc_cover_mat(ii, jj) = sigmoid(cur_w*feature_mat(jj,:,ii).');
        end
    end
end

%% check accuracy
loc_cover_mat

x

% f1 = zeros(n_pt, 1);
% f2 = zeros(n_pt, 1);
%     for i = 1:n_pt
%         % each pt observed sum(loc_cover_mat(i,:)) events
%             f1(i) = x(i) - (provided_label_mat(i,:)*Z_est(2,:).')/(sum(loc_cover_mat(i,:).*Z_est(2,:)) + 1e-4);
%             f2(i) = x(n_pt + i) - (provided_label_mat(i,:)*Z_est(1,:).')/(sum(loc_cover_mat(i,:).*Z_est(1,:)) + 1e-4);
%     end

f = 0;
for j = 1:n_event
    f1 = 0;
    f2 = 0;
    for i = 1:n_pt
        f1 = f1 + provided_label_mat(i,j)*log(loc_cover_mat(i, j)*x(i)) + (1-provided_label_mat(i,j))*log(1-loc_cover_mat(i, j)*x(i));
        f2 = f2 + provided_label_mat(i,j)*log(loc_cover_mat(i, j)*x(n_pt +i)) + (1-provided_label_mat(i,j))*log(1-loc_cover_mat(i, j)*x(n_pt +i));
    end
    f = f + f1*Z_est(2,j) + f2*Z_est(1,j);
end

% Q func is to max; while fminunc() is to min
f = -f;

% add regularization & supervision (to add)
lambda = 80*n_pt;
f = f + lambda*norm(cur_w);


