function loc_cover = loc_cover_mat_by_dis(event_loc, provided_label, thre)

% this function calculates a loc_cover_mat for each participant
% the inputs are 1) event_loc, 2) a pt's provided_label (a row vector)
% it assigns h_ij = 1 if x_ij = 1
% or x_ij = 0 but dis(l_k, L^i) <= thre

% first copy the provided_label since a positive report implies location visit
loc_cover = provided_label;

% the visited locs are those where provide_label == 1
visited_loc = event_loc(logical(provided_label),:);

n = length(provided_label);

for i = 1:n
    if provided_label(i) == 0
        cur_dist = pdist2(event_loc(i,:), visited_loc);
        % min(cur_dist)
        if min(cur_dist) <= thre
            loc_cover(i) = 1;
        end
    end
end
