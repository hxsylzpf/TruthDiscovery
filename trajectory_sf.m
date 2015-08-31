clear;
close all

%% min range - this should be revised upon observing the visiting pattern (where are the locations that pt most frequently visited)

% min_lon = -122.5; max_lon = -122.38;
% min_lat = 37.69; max_lat = 37.81;

% min_lon = -122.44; max_lon = -122.42;
% min_lat = 37.77; max_lat = 37.79;

min_lon = -122.44; max_lon = -122.40;
min_lat = 37.75; max_lat = 37.79;

date_can = [ ...
        2008           5          17
        2008           5          18
        2008           5          19
        2008           5          20
        2008           5          21
        2008           5          22
        2008           5          23
        2008           5          24
        2008           5          25
        2008           5          26
        2008           5          27
        2008           5          28
        2008           5          29
        2008           5          30
        2008           5          31
        2008           6           1
        2008           6           2
        2008           6           3
        2008           6           4
        2008           6           5
        2008           6           6
        2008           6           7
        2008           6           8
        2008           6           9];

% only use 15 days data
date_to_use  = date_can(1:15,:);



% spatial resolution
resol_decimal = 0.001/5; % ~20m
lon_range = min_lon:resol_decimal:max_lon; % col % x-axis
lat_range = min_lat:resol_decimal:max_lat; % row % y-axis

raw_count_mat = sparse(length(lat_range), length(lon_range));
pt_count_mat = sparse(length(lat_range), length(lon_range));

path = 'E:\Datasets\cabspottingdata\data\';

% get all the txt files
list = dir(strcat(path, '*.txt'));

% number of files/pts
n = numel(list);
    
% file format
fmt = '%f %f %f %f';

% maker type
marker_type = {'o', 'x', '>', '*', 'd', '^', 's', 'v'};

% data store
data_store = struct;

% num of participants considered
n_pt = 100;

% the last few pts are chosen to test; while the first few are used to
% build the mobility pattern
n_pt_test = 20;

% rnd seed
seed = 5;

% proportion of hidden
rng(seed);
prop_hidden = 0.4*rand(1, n_pt_test);

% num of event locations sampled
n_sample = 200;

% below how long time to consider a transition is valid
time_interval_limit = 12;

% for all users; all files
for kk = 1:n_pt
    
% for each file, a new local_raw_count_mat is created
local_raw_count_mat = sparse(length(lat_range), length(lon_range));

% initialization - all empty
data_store(kk).lat = [];
data_store(kk).lon = [];
data_store(kk).occu = [];
data_store(kk).time = [];
data_store(kk).loc = [];

        file_name = strcat(path, list(kk).name);
        disp(file_name);  
        % fid = fopen([path, file_name],'r');
        % here, the file_name is already the full path name
        fid = fopen(file_name,'r');
        data = textscan(fid, fmt, 'Delimiter',' ');
        fclose(fid);
        
        % data format
        % 1 Latitude in decimal degree, 2 Longitude in decimal degree, 3 Altitude in feet (-777 if not valid),
        % 4 Date - UNIX time format

        % reverse the trace - time is recorded reversely
        lat = data{1}(end:-1:1);
        lon = data{2}(end:-1:1);
        occu = data{3}(end:-1:1);
        time = data{4}(end:-1:1);
        
        %% round (lon, lat) -> [0.6, 1.4] - round to 1 -> 1 is the center with radius 0.4
        % resolution: 0.0001 approx 11*8m^2 for each grid - very small
        resol = 1/resol_decimal; % 1e3 - convenient for calculation
        lon = round(resol*lon)/resol;
        lat = round(resol*lat)/resol;
        
        % store data as a cell; data for each pt is a cell; each day is concatenated
%         data_store{kk} = [data_store{kk}; [{lat} {lon} {date} {time}]];
        data_store(kk).lat = lat;
        data_store(kk).lon = lon;
        data_store(kk).occu = occu;
        data_store(kk).time = time;
        data_store(kk).loc = [lat lon];
        
        %% plot
%         s = 10;
%         c = linspace(0,1,length(lat));
%         % the last symbol sets the marker type
%         ind = mod(kk, length(marker_type)) + 1;
%         h = scatter(lon,lat,s,c, marker_type{ind});
%         axis([min_lon max_lon min_lat max_lat]);
% %         % annotate
% %         text(lon(1), lat(1), [num2str(kk) ' start'], 'FontSize',14);
% %         text(lon(end), lat(end), [num2str(kk) ' end'], 'FontSize',14);
% %         set(h,'MarkerEdgeColor','b','MarkerFaceColor',[0 .5 .5],'LineWidth',0.6)
%         colorbar
%         hold on

        %% statistics
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % note: if you set as unit8, numerical overflow may happen 
        lon_index = int32((lon - min_lon)*resol + 1);
        lat_index = int32((lat - min_lat)*resol + 1);
        
        % local raw count is for a specific pt % refer to the jth data point
            for j = 1:length(lon_index)
                % only consider the locs that are within the interested range
                if lat_index(j) >= 1 && lon_index(j) >= 1 && lat_index(j) <= length(lat_range) && lon_index(j) <= length(lon_range)
                    local_raw_count_mat(lat_index(j), lon_index(j)) = local_raw_count_mat(lat_index(j), lon_index(j)) + 1;
                end
           end
        
        % update for each user
        raw_count_mat = raw_count_mat + local_raw_count_mat;
        pt_count_mat = pt_count_mat + (local_raw_count_mat > 0);
        
end

%% visualize the two count matrices
%         imagesc(raw_count_mat)
%         colormap summer
        
figure;
imagesc(pt_count_mat);
xlabel('Cell index (x)', 'fontsize', 18);
ylabel('Cell index (y)', 'fontsize', 18);
colormap jet
hcb = colorbar;
set(hcb,'XTickLabel',0:20:100);
set(gca,'YDir','normal', 'fontsize', 18);


% % grids with more than one visit (can be from the same person)
% figure;
% spy(raw_count_mat > 1)
% set(gca,'YDir','normal');

% grids with more than one user visited (from different persons)
figure;
spy(pt_count_mat > 1)
set(gca,'YDir','normal');
title(['Locations where more than 1 pt visited'], 'fontsize', 19);

% imagesc(raw_count_mat)
% colormap gray

% %% construct a sub mat from the original mat
% new_pt_count_mat = pt_count_mat(150:249, 80:179);
% figure;
% imagesc(new_pt_count_mat)
% colormap jet
% colorbar

%% max number of pt visiting an event location
max_p = max(max(pt_count_mat));

% count the number of pt visiting a specific location
% plot the distribution of the number of pt
n_set = 1:max_p;
n_pt_loc = zeros(1, length(n_set));

for k = 1:length(n_set)
    n_pt_loc(k) = sum(sum(pt_count_mat==n_set(k)));
end

n_pt_loc

%% plot pdf and cdf
hf = figure;
% set(hf, 'Position', [200 100 720 420])
hb = bar(n_pt_loc);
xlabel('Number of people visited','fontsize', 18);
ylabel('Number of locations','fontsize', 18);
xlim(gca, [0 100]);
ylim(gca, [0 3000]);

% Save the handle to the axes, hAxes, and get its position.
hAxes = gca;
% set size first and then get the position
set(hAxes, 'XTick', [1 20:20:100], 'YTick', 0:600:3000, 'YGrid', 'on', 'fontsize', 18, 'OuterPosition',[0 0 0.98 1]);
hAxes_pos = get(hAxes,'Position');

% Create a second axes at the same position as the first axes. Store the handle to the second axes, hAxes2. 
hAxes2 = axes('Position',hAxes_pos);

[f,x]= ecdf(pt_count_mat(:));
plot(x, f, 'color', 'g', 'linewidth', 2);
set(hAxes2,'YAxisLocation','right',...
           'Color','none',...
           'XTickLabel',[], 'fontsize', 18);

h1_xlim = get(hAxes,'XLim'); % store x-axis limits of first axes
set(hAxes2,'XLim',h1_xlim); % specify x-axis limits of second axes
ylabel('Cumulative distribution', 'fontsize', 18);

%% sample event locations (index) and plot out these locations
% find index and values of nonzero elements
min_pt = 30; % round(n_pt/4);
linear_idx = find(pt_count_mat >= min_pt);
% [row,col,v] = find(pt_count_mat >= min_pt);
% check the pattern
figure;
spy(pt_count_mat >= min_pt);
set(gca,'YDir','normal');
title(['Locations where more than ' num2str(min_pt) ' pts visited'], 'fontsize', 18);
xlabel('Cell index (x)', 'fontsize', 18);
ylabel('Cell index (y)', 'fontsize', 18);
set(gca, 'fontsize', 16);

% number of those elements
n_total = length(linear_idx);

%% sample event locations

% set seed
% s = rng(seed);
% p = randperm(n_total);
% % location idx
% sample_idx = p(1:n_sample);
linear_idx = linear_idx(randperm(length(linear_idx)));
sample_idx = linear_idx(1:n_sample);
% sample_row_idx = row(p(1:n_sample));
% sample_col_idx = col(p(1:n_sample));

% event locations - note min_lat -> index 1
[I,J] = ind2sub([length(lat_range) length(lon_range)], sample_idx);
sample_lat = min_lat + resol_decimal*(I - 1);
sample_lon = min_lon + resol_decimal*(J - 1);

% to avoid numerical problem - the sampled lat and lon may contain .000001
sample_lat = round(resol*sample_lat)/resol;
sample_lon = round(resol*sample_lon)/resol;

% check pattern
ind_mat = zeros(size(pt_count_mat));
ind_mat(sample_idx) = 1;
figure;
spy(ind_mat);
set(gca,'YDir','normal', 'FontSize', 19);

%% how pt visited these sampled event locations
figure;
plot(sample_lon, sample_lat, 'b*');

event_loc = [sample_lat sample_lon];

% filter out those locations for each pt
n_entry = zeros(n_pt, 1);
n_cpr_loc = zeros(n_pt, 1);

% a structure to store all the reports
report = struct;

for kk = 1:n_pt
    %% retrieve the locations of a pt
    c_loc = [data_store(kk).lat data_store(kk).lon];

    %% find which pt locations are event locations
    % m_idx are indicators 1/0
    % locb are indices wrt to event_loc; 0 means "not a member"
    [m_idx, locb] = ismember(c_loc, event_loc, 'rows');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% plot out both original pt locs, event locs, sampled pt locs to check correctness
%
%     figure;
%     plot(c_loc(:,1), c_loc(:,2), 'bo');
%     hold on
%     plot(event_loc(:,1), event_loc(:,2), 'r*');
%     sampled_c_loc = c_loc(m_idx,:);
%     plot(sampled_c_loc(:, 1), sampled_c_loc(:, 2), 'g+', 'linewidth', 1.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    n_entry(kk) = sum(m_idx);
    
%     report(kk).m_idx = m_idx;
    report(kk).locb = locb(m_idx); % the member index wrt event_loc
    report(kk).pt_loc = c_loc(m_idx,:);
    occu = data_store(kk).occu;
    report(kk).occu = occu(m_idx);
    pt_time = data_store(kk).time;
    report(kk).pt_time = pt_time(m_idx);
    m_time = unixtime2mat(report(kk).pt_time);
    report(kk).pt_time_vec = datevec(m_time);

    %% filter records according to date
    % find the time according to the date_to_use
    time_idx = (zeros(size(report(kk).pt_time,1),1));
    for nn = 1:size(date_to_use, 1)
            time_idx = time_idx | (sum(report(kk).pt_time_vec(:, 1:3) - repmat(date_to_use(nn,:), n_entry(kk), 1), 2) == 0);
    end

    report(kk).locb = report(kk).locb(time_idx); % the member index wrt event_loc
    report(kk).pt_loc = report(kk).pt_loc(time_idx,:);
    report(kk).occu = report(kk).occu(time_idx);
    report(kk).pt_time = report(kk).pt_time(time_idx);
    report(kk).pt_time_vec = report(kk).pt_time_vec(time_idx,:);

    %% trajectory compression
    [report(kk).pt_cpr_loc, idx] = sequence_compress(report(kk).pt_loc);
    % loc index - 1 to n_sample
    report(kk).cpr_locb = report(kk).locb(idx);
    n_cpr_loc(kk) = size(report(kk).pt_cpr_loc, 1);
    report(kk).n_cpr_loc = n_cpr_loc(kk);    
    report(kk).cpr_time = report(kk).pt_time(idx);
    report(kk).cpr_time_vec = report(kk).pt_time_vec(idx,:);
    
    %% find the unique locs and even more compressed expressions
    % compress the results
    [cur_uni, ia, ic] = unique(report(kk).cpr_locb, 'stable');
    n_uni_loc = length(cur_uni);
    % find the time and loc
    report(kk).cpr_uni_time = report(kk).cpr_time(ia);
    report(kk).cpr_uni_time_vec = report(kk).cpr_time_vec(ia,:);
    report(kk).cpr_uni_locb = report(kk).cpr_locb(ia);
    report(kk).cpr_uni_loc = report(kk).pt_cpr_loc(ia,:);
    report(kk).n_cpr_uni_loc = n_uni_loc;
end

save('traj_2008_5_17_21.mat', 'report', 'event_loc');

% save(['traj_' num2str(date_to_use) '.mat'], 'report', 'event_loc');

%% visualize how a pt traveled naturally and visited those event locations
scale = 0;
for kk = 1:2 % n_pt
    tmp1 = report(kk).pt_cpr_loc(:,1);
    tmp2 = report(kk).pt_cpr_loc(:,2);
    figure;
    quiver(tmp1, tmp2, [diff(tmp1); 0], [diff(tmp2); 0], scale);
    hold on
    for aa = 1:length(tmp1)
        text(tmp1(aa), tmp2(aa), num2str(report(kk).cpr_locb(aa)));
        text(tmp1(aa), tmp2(aa)+2*resol_decimal, num2str(report(kk).cpr_time(aa)));
    end
    title(['Participant id' num2str(kk)]);
end

%% build route history for those event locations - from the first few pts; the last few are used for testing

whole_can_record = struct;

for start_id = 1:n_sample
    for end_id = 1:n_sample
        if end_id ~= start_id
            whole_can_record(start_id, end_id).seq = [];
            whole_can_record(start_id, end_id).time = [];
            % consider only the first few pts to build the history
            for kk = 1:n_pt - n_pt_test
                seq = report(kk).cpr_locb;
                time = report(kk).cpr_time;
                % the following function has change the time to relative wrt the 1st record
                [cur_can_seq, cur_can_time] = seq2seg(seq, time, start_id, end_id);
                whole_can_record(start_id, end_id).seq = [whole_can_record(start_id, end_id).seq; cur_can_seq];
                whole_can_record(start_id, end_id).time = [whole_can_record(start_id, end_id).time; cur_can_time];
            end
            % compress the results
            [cur_uni, ia, ic] = unique(whole_can_record(start_id, end_id).seq);
            n_uni = length(cur_uni);
            count = mycount_unique(whole_can_record(start_id, end_id).seq, cur_uni);
            % the unique route sequences
            whole_can_record(start_id, end_id).n_uni = n_uni;
            whole_can_record(start_id, end_id).uni_seq = cur_uni;
            whole_can_record(start_id, end_id).prob_visit = count/sum(count);
            
            % for each unique route; build the histogram of time
            for ii = 1:n_uni
                % logical index
                cur_idx = (ic == ii);
                cur_time = whole_can_record(start_id, end_id).time(cur_idx);
                cur_time = cur_time/3600; % translate into hours
                % build histogram
                [neles, xcenters] = hist(cur_time, 5);
                % prob
                prob_time = neles/sum(neles);
                whole_can_record(start_id, end_id).prob_time(ii) = {num2str(prob_time)};
                whole_can_record(start_id, end_id).time_center(ii) = {num2str(xcenters)};
            end
        
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % check results
% whole_can_record(1,2).seq;
% whole_can_record(1,2).time;
% whole_can_record(1,2).n_uni;
% whole_can_record(1,2).uni_seq;
% whole_can_record(1,2).prob_visit;
% whole_can_record(1,2).prob_time(1);
% whole_can_record(1,2).time_center(1);

% save('loc_visit_hist.mat', 'whole_can_record')

%% given a new user - verify how accurate such estimation is
% again, find how the pt visited event locations
% hide some of the event locations and then verify
% filter out those locations for each pt

% record the est error in loc_cover_mat
err_est_prob = zeros(1, n_pt_test); % prob estimation
err_est_int = zeros(1, n_pt_test); % integer estimation
% % error when assuming all 0 - not necessary - if all are not visited, the
% % model is not complete
% err_0 = zeros(1, n_pt_test);
% error when assuming all 1
err_1 = zeros(1, n_pt_test);
% % error when assuming location range
% err_rad = zeros(1, n_pt_test);
% record the number of revealed visit seq for each participant
n_reveal_cur_visit_seq_record = zeros(1, n_pt_test);


%%
% loc coverage % must be initialized
loc_cover_mat = zeros(n_pt_test, n_sample);

for nn = 1:n_pt_test

        % current pt id
        kk = n_pt - n_pt_test + nn;
        cur_locb = report(kk).cpr_locb;

        [C, ia, ic] = unique(cur_locb);
        
        % build loc_cover_mat
        loc_cover_mat(nn, C) = 1;

        cur_visit_seq = cur_locb(sort(ia));

        length_cur_visit_seq = length(cur_visit_seq);

        n_hidden = round(prop_hidden(nn)*length_cur_visit_seq);

        if n_hidden < 1
            n_hidden = 1;
        end

        % if one of the locs is hidden, at least two should be left
        %%
        if length_cur_visit_seq >= 3
            
           cur_seq_time = report(kk).cpr_time(sort(ia));
           % transform to relative time wrt 1st record
           cur_seq_time = (cur_seq_time - cur_seq_time(1))/3600;
        
            rng(seed);
            cur_p = randperm(length_cur_visit_seq);
            % copy the original seq
            reveal_cur_visit_seq = cur_visit_seq;
            reveal_cur_seq_time = cur_seq_time;
            % make the hidden ones deleted
            reveal_cur_visit_seq(cur_p(1:n_hidden)) = [];
            reveal_cur_seq_time(cur_p(1:n_hidden)) = [];

            n_reveal_cur_visit_seq = length(reveal_cur_visit_seq);

            n_reveal_cur_visit_seq_record(nn) = n_reveal_cur_visit_seq;
            
            % predict the hidden ones - since it is a sequence, only consider the latter one as the destination
            outer_record_loc = [];
            outer_record_prob = [];

            %%
            for aa = 1:n_reveal_cur_visit_seq-1
                %%
                for bb = aa+1:n_reveal_cur_visit_seq
                    %%
                    if bb ~= aa
                        cur_time = reveal_cur_seq_time(bb) - reveal_cur_seq_time(aa);
                        % id of the starting location
                        start_id = reveal_cur_visit_seq(aa);
                        end_id = reveal_cur_visit_seq(bb);

                        n_route = whole_can_record(start_id,end_id).n_uni;
                        %%
                        if n_route > 0

                        % find the probability of taking a route
                            cur_prob_visit = whole_can_record(start_id,end_id).prob_visit;
                            % assume a little basis prob
                            cur_prob_visit_time = 1e-2*ones(size(cur_prob_visit));
                        % find the probability of using current time
                            for i = 1:n_route
                                % use {} here since the stored info is in cell format
                                cur_time_center = str2num(whole_can_record(start_id,end_id).time_center{i});
                                cur_prob_time = str2num(whole_can_record(start_id,end_id).prob_time{i});

                                cur_time_diff = abs(cur_time_center - cur_time);
                                [min_diff, center_idx] = min(cur_time_diff);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % ignore time that is larger than 24 hours
                                if cur_time_center(center_idx) <= time_interval_limit
                                    cur_prob_visit_time(i) = cur_prob_time(center_idx);
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            end
                        % final probability of taking a route within current time
                            cur_final_prob = cur_prob_visit.*cur_prob_visit_time;
                        % normalize
                            cur_final_prob = cur_final_prob/sum(cur_final_prob);

                        % merge the probabilities
                            idx = (cur_final_prob > 0);
                            retain_route = whole_can_record(start_id,end_id).uni_seq(idx);
                            retain_final_prob = cur_final_prob(idx);

                         % routes that remains; with a prob larger than 0
                            n_retain_route = sum(idx);
                            record_uni_route = [];
                            record_prob = [];
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            %% only continue when n_retain_route is at least 1
                            if n_retain_route > 0
                                for cc = 1:n_retain_route
                                    cur_route = str2num(retain_route{cc});
                                    cur_uni_route = setdiff(unique(cur_route), reveal_cur_visit_seq);
                                    record_uni_route = [record_uni_route cur_uni_route];
                                    record_prob = [record_prob repmat(retain_final_prob(cc), 1, length(cur_uni_route))];
                                end
                                [merged_uni_loc, ~, ic] = unique(record_uni_route);
                                 
                                merged_prob = zeros(size(merged_uni_loc));
                                
                                for dd = 1:length(merged_uni_loc)
                                    merged_prob(dd) = sum(record_prob(record_uni_route == merged_uni_loc(dd)));
                                end
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % outer record for each start and end location combination
                                outer_record_loc = [outer_record_loc merged_uni_loc];
                                outer_record_prob = [outer_record_prob merged_prob];
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                            end
                            %
                        end

                    end

                end
            end
            % finally take the max probability as the probability of visiting a location
                     [final_uni_loc, ia, ic] = unique(outer_record_loc);
                     % must initiate it for every new computation
                     final_prob = [];
                     for ee = 1:length(final_uni_loc)
                            final_prob(ee) = max(outer_record_prob(outer_record_loc == final_uni_loc(ee)));
                     end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%
        % check results
        cur_visit_seq
        reveal_cur_visit_seq
        final_uni_loc
        final_prob

        %% compute estimation error
        % estimation is made based on those unobserved
        est_loc = setdiff(1:n_sample, reveal_cur_visit_seq);
        n_est_loc = length(est_loc);
        true_visit_vec = zeros(1, n_est_loc);
        est_visit_vec = zeros(1, n_est_loc);

        % construct the ground truth and est vector
        for ff = 1:n_est_loc
            cur_est_loc = est_loc(ff);
            % cur_visit_seq is the true visit seq
            idx1 = find(cur_visit_seq == cur_est_loc);
            % if found
            if sum(idx1) > 0
                true_visit_vec(ff) = 1;
            end

            % final_uni_loc is the est visit loc
            idx2 = find(final_uni_loc == cur_est_loc);
            if sum(idx2) > 0
                est_visit_vec(ff) = final_prob(idx2);
            end

        end

        % calculate the error % use Ave or RMSE
        err_est_prob(nn) = mean(abs(true_visit_vec - est_visit_vec));
        % round the prob est to integer
        err_est_int(nn) = mean(abs(true_visit_vec - round(est_visit_vec)));
%         err_0(nn) = mean(abs(true_visit_vec - zeros(size(true_visit_vec))))
        err_1(nn) = mean(abs(true_visit_vec - ones(size(true_visit_vec))));
        
        % err(nn) = sqrt(mean((true_visit_vec - est_visit_vec).^2))
end

% final average error
mean_err = [1-mean(err_1) mean(err_est_int) mean(err_est_prob) ]
std_err = [std(err_1) std(err_est_int) std(err_est_prob)]

x_val = n_reveal_cur_visit_seq_record/n_sample;

plot(x_val, err_1, 'rv', x_val, err_est_int, 'gx', x_val, err_est_prob, 'bo', 'linewidth', 1.5);
xlabel('Proportion of reported locations', 'fontsize', 18);
ylabel('Location visit estimation accuracy', 'fontsize', 18);
set(gca, 'fontsize', 18);

% loc_cover stat
n_pt_per_loc = sum(loc_cover_mat, 1)
n_loc_per_pt = sum(loc_cover_mat, 2).'


