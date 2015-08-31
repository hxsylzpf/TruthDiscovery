function [f, grad] = func_obj_grad_abw_event(x, n_pt, n_event, n_feature, provided_label_mat, feature_mat, Z_est)

% feature_mat = zeros(n_event, n_feature, n_pt)

% x is constructed as [a b w d]
% a - x(1) to x(n_pt)
% b - x(n_pt + 1) to x(2*n_pt)
% w - x(2*n_pt + 1) to x(2*n_pt + n_event*n_feature)
% w_j - x(2*n_pt + (i-1)*n_feature + 1) to x(2*n_pt + i*n_feature)
% d - x(2*n_pt + n_event*n_feature + 1)

% construct the loc_cover_mat - first copy ones from provided_label_mat
loc_cover_mat = provided_label_mat;

% then construct prob est according to logistic regression
cur_w = zeros(n_feature, n_event);
% time_thre = 8; % hours
for jj = 1:n_event
    cur_w(:,jj) = x(2*n_pt + (jj-1)*n_feature + 1 : 2*n_pt + jj*n_feature).';
     for ii = 1:n_pt
        if loc_cover_mat(ii, jj) == 0 % && feature_mat(jj,2,ii) < time_thre % not a revealed loc & time (2nd feature) spent is less than 4
            loc_cover_mat(ii, jj) = sigmoid(feature_mat(jj,:,ii)*cur_w(:,jj));
%         elseif loc_cover_mat(ii, jj) == 0 && feature_mat(jj,2,ii) >= time_thre % not a revealed loc & time spent is larger than 4
%             loc_cover_mat(ii, jj) = 0.5;
        end
    end
end

% check accuracy
loc_cover_mat(1:10, 1:5)

% x

%% obj function
f = 0;
for j = 1:n_event
    f1 = 0;
    f2 = 0;
    for i = 1:n_pt
        f1 = f1 + provided_label_mat(i,j)*mylog(loc_cover_mat(i,j)*x(i)) + (1-provided_label_mat(i,j))*mylog(1-loc_cover_mat(i,j)*x(i));
        f2 = f2 + provided_label_mat(i,j)*mylog(loc_cover_mat(i,j)*x(n_pt +i)) + (1-provided_label_mat(i,j))*log(1-loc_cover_mat(i,j)*x(n_pt +i));
    end
    f1 = f1 + mylog(x(2*n_pt + n_event*n_feature + 1));
    f2 = f2 + mylog(1 - x(2*n_pt + n_event*n_feature + 1));
    f = f + f1*Z_est(2,j) + f2*Z_est(1,j);
end

% Q func is to max; while fminunc() is to min
f = -f;

% add regularization - minimize the regularization
lambda = n_event;
for jj = 1:n_event
    f = f + lambda*norm(cur_w(:,jj), 1);   
end

%% gradient
if nargout > 1
    grad_a = zeros(n_pt, 1);
    grad_b = zeros(n_pt, 1);
    grad_w = zeros(n_feature, n_event); % each event has a length-n_feature w_j
    grad_d = 0;
    
    for i = 1:n_pt
        for j = 1:n_event
            grad_a(i) = grad_a(i) + Z_est(2,j)*(provided_label_mat(i,j) - loc_cover_mat(i, j)*x(i))/(x(i)*(1 - loc_cover_mat(i, j)*x(i)));
            grad_b(i) = grad_b(i) + Z_est(1,j)*(provided_label_mat(i,j) - loc_cover_mat(i, j)*x(n_pt+i))/(x(n_pt+i)*(1 - loc_cover_mat(i, j)*x(n_pt+i)));    
        end    
    end

    for j = 1:n_event
        grad_d = grad_d + Z_est(2,j)/x(2*n_pt + n_event*n_feature + 1) - Z_est(1,j)/(1-x(2*n_pt + n_event*n_feature + 1));
        for i = 1:n_pt
             if loc_cover_mat(i,j) ~= 1
             grad_w(:,j) = grad_w(:,j) + Z_est(2,j)*(provided_label_mat(i,j) - loc_cover_mat(i, j)*x(i))*(1 - loc_cover_mat(i, j))/(1 - loc_cover_mat(i, j)*x(i))*feature_mat(j,:,i).' ...
                                   + Z_est(1,j)*(provided_label_mat(i,j) - loc_cover_mat(i, j)*x(n_pt+i))*(1 - loc_cover_mat(i, j))/(1 - loc_cover_mat(i, j)*x(n_pt+i))*feature_mat(j,:,i).';
            end     
        end
    end
    
    grad = [grad_a; grad_b; grad_w(:); grad_d];
    grad = -grad;
end

