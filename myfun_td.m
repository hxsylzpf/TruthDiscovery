function F = myfun_td(x, n_pt, n_event, provided_label_mat, loc_cover_mat, z_post)

% x is constructed as [alpha beta0]
% alpha - x(1) to x(n_participant)
% beta0 - x(n_participant + 1) to x(n_participant + n_event)

f1 = zeros(n_pt, 1);
for i = 1:n_pt
    for j = 1:n_event
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update considering only 1s but not 0s
        if loc_cover_mat(i,j) == 1
%         %%%%%%%%%%%%%%%%%%%%%%%%%%
        f1(i) = f1(i) + exp(x(n_pt + j))*(z_post(2,j)*provided_label_mat(i,j) + z_post(1,j)*(1 - provided_label_mat(i,j)) - 1/(1+exp(-x(i)*exp(x(n_pt + j)))));
        end
    end
end

f2 = zeros(n_event, 1);
for j = 1:n_event
    for i = 1:n_pt
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % update considering only 1s but not 0s
        if loc_cover_mat(i,j) == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        f2(j) = f2(j) + x(i)*(z_post(2,j)*provided_label_mat(i,j) + z_post(1,j)*(1 - provided_label_mat(i,j)) - 1/(1+exp(-x(i)*exp(x(n_pt + j)))));
        end
    end
end

F = [f1; f2];
