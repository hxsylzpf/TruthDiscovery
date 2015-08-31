clear;
close all

%% min range - this should be revised upon observing the visiting pattern (where are the locations that pt most frequently visited)

% min_lon = -122.5; max_lon = -122.38;
% min_lat = 37.69; max_lat = 37.81;

min_lon = -122.44; max_lon = -122.40;
min_lat = 37.75; max_lat = 37.79; % area1:75-77, area2:77-79

% min_lon = -122.44; max_lon = -122.40;
% min_lat = 37.77; max_lat = 37.79;

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

% % only use 15 days data
% date_to_use  = date_can(1:15,:);

% spatial resolution
resol_decimal = 0.001/2; % ~50m % 0.00001; %
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
n_pt = 200;
n_uni_day_record = zeros(n_pt,1);

% a pt should have 15 day's data
min_n_day = 15;

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
        
        % find the number of days in each pt's data
        
        unix_time = unixtime(time);
        uni_day = unique(unix_time(:,[2 3]), 'rows');
        n_uni_day = size(uni_day,1);
        n_uni_day_record(kk) = n_uni_day;
        
        
        %% round (lon, lat) -> [0.6, 1.4] - round to 1 -> 1 is the center with radius 0.4
        % resolution: 0.0001 approx 11*8m^2 for each grid - very small
        resol = 1/resol_decimal; % 1e3 - convenient for calculation
        lat = round(resol*lat)/resol;
        lon = round(resol*lon)/resol;
                
        % get the idx of locs that fall inside the area of interest
        idx_in_range = (lat >= min_lat) & (lat <= max_lat) & (lon >= min_lon) & (lon <= max_lon);

        % record only when the pt has some data inside the area of interest
        if sum(idx_in_range) > 0  
            lat = lat(idx_in_range);
            lon = lon(idx_in_range);
            occu = occu(idx_in_range);
            time = time(idx_in_range);

%           store data as struct; data for each pt is a struct; each day is concatenated
            data_store(kk).lat = lat;
            data_store(kk).lon = lon;
            data_store(kk).occu = occu;
            data_store(kk).time = time;
            data_store(kk).loc = [lat lon];
        else
            data_store(kk).lat = [];
            data_store(kk).lon = [];
            data_store(kk).occu = [];
            data_store(kk).time = [];
            data_store(kk).loc = [];    
        end
        
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

n_uni_day_record

qualified_pt_idx = (n_uni_day_record >= min_n_day);

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
title('Locations where more than 1 pt visited', 'fontsize', 19);

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
xlim(gca, [0 200]);
ylim(gca, [0 4000]);

% Save the handle to the axes, hAxes, and get its position.
hAxes = gca;
% set size first and then get the position
set(hAxes, 'XTick', [1 25:25:200], 'YTick', 0:400:4000, 'XGrid', 'on', 'YGrid', 'on', 'fontsize', 18, 'OuterPosition',[0 0 0.98 1]);
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

%% process each pt's location trace
    %% find the unique locs

    time_traffic_light = [20 120];
    dist_traffic_light = [0 70];
    
    n_total_report = 0;
    all_report_loc = [];
    
%% time and distance based processing
% for kk = 1:n_pt
%     
%     % compress the results    
%     [data_store(kk).uni_loc, ia] = sequence_compress(data_store(kk).loc);
%     
% %     figure;
% %     plot(data_store(kk).loc(:,1), data_store(kk).loc(:,2), 'b');
% %     figure;
% %     plot(data_store(kk).uni_loc(:,1), data_store(kk).uni_loc(:,2), 'g');
% 
%     n_uni_loc = length(ia);
%     % get the unique time
%     data_store(kk).uni_time = data_store(kk).time(ia);
%     
%     n_diff = n_uni_loc - 1;
%     data_store(kk).time_diff = zeros(n_diff, 1);
%     data_store(kk).dist = zeros(n_diff, 1);
% 
%     for i = 1:n_diff
%         % time difference % in second
%         data_store(kk).time_diff(i) = data_store(kk).uni_time(i+1) - data_store(kk).uni_time(i);
%         [d1km, ~]=lldistkm(data_store(kk).uni_loc(i,:),data_store(kk).uni_loc(i+1,:));
%         % distance % in meter
%         data_store(kk).dist(i) = d1km*1e3;
%     end
% 
%     %% process reports
%     % filter based on time diff and dist
%     report_idx = (data_store(kk).time_diff >= time_traffic_light(1)) & (data_store(kk).time_diff <= time_traffic_light(2)) ...
%         & (data_store(kk).dist <= dist_traffic_light(2));
%     data_store(kk).n_report = sum(report_idx)
%     n_total_report = n_total_report + data_store(kk).n_report;
%     
%     %% store the time when reporting
%     data_store(kk).report_time = data_store(kk).uni_time([report_idx; false]);   
%     % record report location
%     data_store(kk).report_loc = data_store(kk).uni_loc([report_idx;false],:);
%     all_report_loc = [all_report_loc; data_store(kk).report_loc];
% end

%% velocity and time based processing
for kk = 1:n_pt
    
    for i = 1:length(data_store(kk).time)-1
        % time difference % in second
        data_store(kk).time_diff(i) = data_store(kk).time(i+1) - data_store(kk).time(i);
        [d1km, ~]=lldistkm(data_store(kk).loc(i,:),data_store(kk).loc(i+1,:));
        % distance % in meter
        data_store(kk).dist(i) = d1km*1e3;
        data_store(kk).vel(i) = data_store(kk).dist(i)/data_store(kk).time_diff(i);
    end
    
    % % [idx22, dur_signal] = cum_time_for_segment(idx11, data_store(kk).time(1:end-1));

    %% detection method for uninterpolated data
    % find point where vel < 1m/s (such that the pt is waiting)
    report_idx = (data_store(kk).vel < 3) & (data_store(kk).time_diff >= 20) & (data_store(kk).time_diff <= 180);
    data_store(kk).n_report = sum(report_idx)
    n_total_report = n_total_report + data_store(kk).n_report;

    %% store the time when reporting % record both the start and end time/loc
    %% detection method for uninterpolated data
    % interpolate the data between each reported_idx locs and record
    data_store(kk).report_time = [data_store(kk).time([report_idx false]); data_store(kk).time([false report_idx])];    
    % % record report location
    % data_store(kk).report_loc = [data_store(kk).loc([report_idx false],:); data_store(kk).loc([false report_idx],:)];
    data_store(kk).report_loc = data_store(kk).loc([false report_idx],:);
    all_report_loc = [all_report_loc; data_store(kk).report_loc];
end

n_total_report

uni_report_loc = unique(all_report_loc, 'rows');
n_uni_report_loc = size(uni_report_loc,1)

%% plot several pt's data
my_color = {'bo','gx','r+','c*','mv','y>'};
delta = 5;
figure;
hold on
for ii = 1:6
    plot(data_store(delta*ii).report_loc(:,1), data_store(delta*ii).report_loc(:,2), my_color{ii});
end
hold off
xlim([min_lat max_lat]);
ylim([min_lon max_lon]);

    % if the data set contain samples larger than the value of the last bin,
    % all of them will be put in the last bin
    xcen_time = [2 5 10:10:200];
    nele_time = hist(data_store(kk).time_diff, xcen_time);
    figure;
    bar(xcen_time, nele_time);

    xcen_dist = 10:10:500;

    nele_dist = hist(data_store(kk).dist, xcen_dist);
    figure;
    bar(xcen_dist, nele_dist);

% %% clustering the reported locations
% % resol_decimal2 = 0.001/5;
% % resol = 1/resol_decimal2; % 1e3 - convenient for calculation
% % cpr_lat = round(resol*all_report_loc(:,1))/resol;
% % cpr_lon = round(resol*all_report_loc(:,2))/resol;
% 
% % simplest way is to find unique
% uni_report_loc = unique(all_report_loc,'rows');
% n_uni_report_loc = size(uni_report_loc,1);
% 
% figure;
% plot(uni_report_loc(:,1), uni_report_loc(:,2), 'bo');
% xlim([min_lat max_lat]);
% ylim([min_lon max_lon]);

%% clustering the reported locations
n_dim = 2;
Y = pdist(uni_report_loc, 'cityblock');
Z = linkage(Y, 'average'); % 'average', 'complete'
dendrogram(Z);
c = cophenet(Z,Y)

%%
cutoff = 0.001/3;

T = cluster(Z,'cutoff',cutoff,'criterion','distance');

n_cluster = length(unique(T));

% calculate cluster centers
cluster_center = zeros(n_cluster, n_dim);
n_ele_per_cluster = zeros(n_cluster, 1);
for i = 1:n_cluster
    cur_idx = (T == i);
    n_ele_per_cluster(i) = sum(cur_idx);
    cluster_center(i,:) = median(uni_report_loc(cur_idx,:), 1);
end

% idx_valid_cluster = (n_ele_per_cluster >= 1);
% cluster_center = cluster_center(idx_valid_cluster,:);
n_cluster = size(cluster_center, 1);

% visualize all the data and cluster centers
% my_color = {'bo', 'rx', 'c*', 'm>', 'bd', 'y+', 'rx', 'gv', 'm*', 'bs', 'c+', 'cv', 'go', 'yx', 'g*', 'g>', 'yd', 'y+', 'bx', 'rv', 'm*', 'ms', 'c+', 'yv', ...
%                         'bo', 'rx', 'c*', 'm>', 'bd', 'y+', 'rx', 'gv', 'm*', 'bs', 'c+', 'cv', 'go', 'yx', 'g*', 'g>', 'yd', 'y+', 'bx', 'rv', 'm*', 'ms', 'c+', 'yv'};
% my_color = repmat(my_color, 1, 10);
% scrsz = get(0,'ScreenSize');
% % scrsz(3) - width; scrsz(4) - height
% figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/1.5 scrsz(3)/4]);
% h = 0;
% hold on
% for i = 1:n_cluster
%         cur_idx = (T == i);
%         plot(uni_report_loc(cur_idx,1), uni_report_loc(cur_idx,2), my_color{i}, 'linewidth', 1.5);
%         plot(cluster_center(i,1), cluster_center(i,2), my_color{i}, 'markersize', 14);
%         text(cluster_center(i,1) + 0.02, cluster_center(i,2), num2str(i), 'fontsize', 20);
% end
% hold off
% 
% xlabel('Lat', 'fontsize', 20);
% ylabel('Lon', 'fontsize', 20);
% % xlim([1.2 4]); ylim([0.7 2.5])
% % set(gca, 'fontsize', 20, 'XTick', [1.2 1.5:0.5:4]);

%% randomly pick out 100 cluster centers
n_cluster_to_use = 1000; % 100;
seed = 1;
% rng(seed);
idx_perm = randperm(n_cluster);

cluster_center_to_use = cluster_center(idx_perm(1:n_cluster_to_use),:);
figure;
plot(cluster_center_to_use(:,1), cluster_center_to_use(:,2), 'bo');
grid on;

% save these cluster centers to file

%% read in the locations w ground truth

% % read in event label
% % fileID = fopen('report_loc_area1.csv');
% % C = textscan(fileID, '%f%f%f%*s%f%f%f%*s%f', 'delimiter', ',');
% fileID = fopen('report_loc_area1_radius50_20141014.csv');
% C = textscan(fileID, '%f%f%f%*f%*f%*f%*f', 'delimiter', ',');
% fclose(fileID);
% anno_loc = [C{1} C{2}];
% event_label = C{3};
% % lam_g1 = C{4};
% % lam_g0 = n_pt - lam_g1;
% % loc_popu = lam_g1/n_pt;
% % 
% % n_anno_loc = size(anno_loc,1);
% % 
% % filename = 'report_loc_SF_area1_20140915_anno.xlsx';
% % sheet = 1;
% % xlRange = 'A2:G151';
% % 
% % data = xlsread(filename, sheet, xlRange);
% % 
% % anno_loc = data(:,1:2);
% % % n_pt_per_event = data(:,3);
% % % n_visit_per_event = data(:,4);
% % % report_ratio = data(:,5);
% % % est_label = data(:,6);
% % event_label = data(:,7);

%% find who reported which loc - provided_label_mat and loc_cover_mat

% % if use the cluster center directly
anno_loc = cluster_center_to_use;
n_anno_loc = size(anno_loc,1);

provided_label_mat = zeros(n_pt, n_anno_loc);
loc_cover_mat = zeros(n_pt, n_anno_loc);
% thre for considering as close enough
thre = 0.001/3; % 50m

% record the number of report and pass by
provided_label_mat_num = zeros(n_pt, n_anno_loc);
loc_cover_mat_num = zeros(n_pt, n_anno_loc);

for kk = 1:n_pt
%     % whether the aggregated unique reported locs can be found in the reported locs
%     [Lia,Locb] = ismember(anno_loc, data_store(kk).report_loc,'rows');
%     provided_label_mat(kk, Lia) = 1;
    % whether the reported locations can be found in the original locations
%     [Lia1, Locb1] = ismember(anno_loc, data_store(kk).loc,'rows');
%     loc_cover_mat(kk,Lia1) = 1;

%% do not take into account multiple reports at the same location
    if ~isempty(data_store(kk).report_loc)
        % compute distance to each annotated loc
        [D,~] = pdist2(data_store(kk).report_loc, anno_loc,'cityblock','Smallest',1);
        provided_label_mat(kk, :) = (D <= thre);
        % compute distance to each annotated loc
        [D,~] = pdist2(data_store(kk).loc, anno_loc,'cityblock','Smallest',1);
        loc_cover_mat(kk, :) = (D <= thre);
    end

%% consider multiple reports at the same location
% if ~isempty(data_store(kk).report_loc)
%     % for each loc with gt, find how many corresponding reported locs
%     for i = 1:n_anno_loc
%         temp_dist1 = pdist2(data_store(kk).report_loc, anno_loc(i,:));
%         provided_label_mat_num(kk,i) = sum(temp_dist1 <= thre);
%         
%         % calculate loc visit based on days
%         temp_dist2 = pdist2(data_store(kk).loc, anno_loc(i,:));
%         % % remove the count due to the repetition of fixed-interval sampling
%         % % e.g., [1 1 1 2 2 2 3 3 3]
%         % temp_dist2 = sequence_compress(temp_dist2);
%         loc_cover_mat_num(kk,i) = sum(temp_dist2 <= thre);
%     end    
% end

end

%% process num mat to get binary mat
% loc_cover_mat_new = (loc_cover_mat_num >= 1) | (provided_label_mat > 0);
% 
% figure;
% spy(loc_cover_mat)
% figure;
% spy(loc_cover_mat_new)
% 
% loc_cover_mat = loc_cover_mat_new;

%% check correctness
diff = loc_cover_mat - provided_label_mat;

if sum(diff(:) < 0) == 0
    disp('Correct!');
else
    disp('Error!');
end

n_pt_per_event = sum(provided_label_mat).'

max_n_pt = max(n_pt_per_event)

n_visit_per_event = sum(loc_cover_mat).'

report_ratio = n_pt_per_event./n_visit_per_event

est_label = (report_ratio >= 0.5);

sum(est_label)

%% 
sum(event_label)

err = sum(est_label ~= event_label)

idx_correct = (est_label == event_label);

%% save data 
% % lat, lon, n_report, n_visit, ratio
% M = [cluster_center_to_use n_pt_per_event n_visit_per_event report_ratio est_label];
% xlswrite('cluster_cen_area1_visit2.xlsx', M);

%% record the report loc and time for those final event locs
for kk = 1:n_pt
    % whether report
    [Lia,Locb] = ismember(data_store(kk).report_loc, anno_loc,'rows');
    data_store(kk).report_loc_filter = data_store(kk).report_loc(Lia,:);
    data_store(kk).report_time_filter = data_store(kk).report_time(Lia);
end

%% filter based on mimimum number of pts and evts
min_n_pt = 10;
min_n_evt = 10;
[provided_label_mat_filter, loc_cover_mat_filter, report_loc_filter, idx_evt_w_enough_pt] = ...
    n_pt_n_evt_filter(min_n_pt, min_n_evt, provided_label_mat, loc_cover_mat, anno_loc);

[n_pt,n_event] = size(provided_label_mat_filter);

event_label = event_label(idx_evt_w_enough_pt).';

lam_g1 = sum(loc_cover_mat_filter).';
lam_g0 = n_pt - lam_g1;

loc_popu = lam_g1/n_pt;
% lam_g1 = lam_g1(idx_evt_w_enough_pt);
% lam_g0 = lam_g0(idx_evt_w_enough_pt);

%% spy the two mats
figure;
spy(provided_label_mat_filter)
title('Provided_label');

figure;
spy(loc_cover_mat_filter)
title('Loc coverage');

%% location cover ratio of the participants
loc_cover_ratio = sum(loc_cover_mat_filter,2)/n_event;

figure;
stem(sort(loc_cover_ratio));
title('Location coverage per pt');

%% location popularity
loc_popularity = sum(loc_cover_mat_filter,1)/n_pt;
figure;
stem(sort(loc_popularity));
title('Location popularity per evt');

%% response ratio
res_ratio = sum(provided_label_mat_filter,2)./sum(loc_cover_mat_filter,2);
figure;
stem(sort(res_ratio));
title('Response ratio');

%% final locations
figure;
plot(report_loc_filter(:,1), report_loc_filter(:,2), 'bo');
grid on;

figure;
n_pt_per_event_filter = sum(provided_label_mat_filter,1);
nele = hist(n_pt_per_event_filter,min_n_pt:max_n_pt);
bar(min_n_pt:max_n_pt, nele);

% save to file
save('data_traffic_light_SF_pt200_evt160_area1_radius50_mdf_20141014.mat','provided_label_mat_filter', 'loc_cover_mat_filter', 'event_label', ...
    'data_store', 'n_pt', 'n_event', 'report_loc_filter', 'lam_g1', 'lam_g0', 'loc_popu');
