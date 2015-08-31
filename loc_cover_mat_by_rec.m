function loc_cover = loc_cover_mat_by_rec(event_loc, provided_label)

% this function calculates a loc_cover_mat for each participant
% it assumes a rectangle where the pt visited all the locs inside it
% but did not visit all the locs outside it

% first copy the provided_label since a positive report implies location visit
loc_cover = provided_label;

n = length(provided_label);

revealed_loc = event_loc(logical(provided_label),:);

% find boundaries
min_lat = min(revealed_loc(:,1));
max_lat = max(revealed_loc(:,1));
min_lon = min(revealed_loc(:,2));
max_lon = max(revealed_loc(:,2));

for i = 1:n
    if (provided_label(i) == 0) && (event_loc(i,1) >= min_lat) && (event_loc(i,1) <= max_lat) && ...
            (event_loc(i,2) >= min_lon) && (event_loc(i,2) <= max_lon)
        loc_cover(i) = 1;
    end
end

%% check correctness
% figure;
% plot(event_loc(:,1), event_loc(:,2), 'k.', ... % all event locs
%      event_loc(logical(loc_cover), 1), event_loc(logical(loc_cover), 2), 'go', ... % ested coverage
%      event_loc(logical(provided_label), 1), event_loc(logical(provided_label), 2), 'b*', ... % 
%      [min_lat max_lat], [min_lon max_lon], 'rd'); % boundary points

