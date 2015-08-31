clear;
close all

%% min range - this should be revised upon observing the visiting pattern (where are the locations that pt most frequently visited)
% % original range
% min_lon = 116.2; max_lon = 116.45;
% min_lat = 39.8; max_lat = 40.05;

% % revised range
min_lon = 116.28; max_lon = 116.38;
min_lat = 39.95; max_lat = 40.05;

resol_decimal = 0.001;
lon_range = min_lon:resol_decimal:max_lon; % col % x-axis
lat_range = min_lat:resol_decimal:max_lat; % row % y-axis

% not very large - 251*251
% # of rows (rows - lat) * # of columns (cols - lon)
raw_count_mat = sparse(length(lat_range), length(lon_range));
pt_count_mat = sparse(length(lat_range), length(lon_range));

username_set = 0:14;
n_pt = length(username_set);

% data store
data_store = struct;

%% for all the users
for kk = 1:n_pt
    
% username = '001';
% formating username
username = sprintf('%03d', username_set(kk));
path = ['E:\Datasets\Geolife Trajectories 1.3\Data_subset\', username, '\Trajectory\'];

%% for each user - for all the trajectories
% get all the plt files
list = dir(strcat(path, '*.plt'));

% number of files
n = numel(list);
    
% file format
fmt = '%f %f %*f %f %s %s %s';
    
% maker type
marker_type = {'o', 'x', '>', '*', 'd', '^', 's', 'v'};

% for each file % or for a subset of files

% for each user, a new local_raw_count_mat is created
local_raw_count_mat = sparse(length(lat_range), length(lon_range));

% initialization - all empty
data_store(kk).lat = [];
data_store(kk).lon = [];
data_store(kk).date = '';
data_store(kk).time = '';

for i = 1:n
    
        file_name = strcat(path, list(i).name);
        disp(file_name);  
        % fid = fopen([path, file_name],'r');
        % here, the file_name is already the full path name
        fid = fopen(file_name,'r');
        data = textscan(fid, fmt, 'Delimiter',',', 'Headerlines',6);
        
        % data format
        % 1 Latitude in decimal degree, 2 Longitude in decimal degree, 3 Altitude in feet (-777 if not valid),
        % 4 Date - number of days (with fractional pt) that have passed since 12/30/1899, 5 Date as a string, 6 Time as a string

        lat = data{1};
        lon = data{2};
        date = data{5}; % if in the same file, the date is generally the same day
        time = data{6};
        
        %% round (lon, lat) -> [0.6, 1.4] - round to 1 -> 1 is the center with radius 0.4
        % resolution: 0.0001 approx 11*8m^2 for each grid - very small
        resol = 1/resol_decimal; % 1e3 - convenient for calculation
        lon = round(resol*lon)/resol;
        lat = round(resol*lat)/resol;
        
        % store data as a cell; data for each pt is a cell; each day is concatenated
%         data_store{kk} = [data_store{kk}; [{lat} {lon} {date} {time}]];
        data_store(kk).lat = [data_store(kk).lat; lat];
        data_store(kk).lon = [data_store(kk).lon; lon];
        data_store(kk).date = [data_store(kk).date; date];
        data_store(kk).time = [data_store(kk).time; time];
        
        %% plot
%         s = 10;
%         c = linspace(0,1,length(lat));
%         % the last symbol sets the marker type
%         ind = mod(i, length(marker_type)) + 1;
%         h = scatter(lon,lat,s,c, marker_type{ind});
%         % annotate the 
% %         text(lon(1), lat(1), [num2str(i) ' start'], 'FontSize',14);
% %         text(lon(end), lat(end), [num2str(i) ' end'], 'FontSize',14);
%         % set(h,'MarkerEdgeColor','b','MarkerFaceColor',[0 .5 .5],'LineWidth',0.6)
%         % colorbar
%         hold on

        %% statistics
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % note: if you set as unit8, numerical overflow may happen 
        lon_index = uint32((lon - min_lon)*resol + 1);
        lat_index = uint32((lat - min_lat)*resol + 1);
        
        % local raw count is for a specific pt
        for j = 1:length(lon_index)
            if lat_index(j) >= 1 && lon_index(j) >= 1 && lat_index(j) <= length(lat_range) && lon_index(j) <= length(lon_range)
                local_raw_count_mat(lat_index(j), lon_index(j)) = local_raw_count_mat(lat_index(j), lon_index(j)) + 1;
            end
        end     
        
end
% update for each user
        raw_count_mat = raw_count_mat + local_raw_count_mat;
        pt_count_mat = pt_count_mat + (local_raw_count_mat > 0);
end
% grid on
% xlim([116.2 116.45]);
% ylim([39.8 40.05]);
% 
% xlabel('Longitude', 'fontsize', 20);
% ylabel('Latitude', 'fontsize', 20);
% set(gca, 'FontSize', 19);

%% visualize the two count matrices
%         imagesc(raw_count_mat)
%         colormap summer
        
figure;
imagesc(pt_count_mat)
colormap gray % jet
% colorbar
set(gca,'YDir','normal');

% grids with more than one visit (can be from the same person)
figure;
spy(raw_count_mat > 1)
set(gca,'YDir','normal');

% grids with more than one user visited (from different persons)
figure;
spy(pt_count_mat > 1)
set(gca,'YDir','normal');

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

figure;
bar(n_pt_loc)
xlabel('Number of participants visited','fontsize', 20);
ylabel('Count of locations','fontsize', 20);
set(gca, 'FontSize', 19);

%% sample event locations (index) and plot out these locations
% find index and values of nonzero elements
min_pt = 3;
[row,col,v] = find(pt_count_mat >= min_pt);
% check the pattern
figure;
spy(pt_count_mat >= min_pt);
set(gca,'YDir','normal');
title(['Locations where more than ' num2str(min_pt) ' pt visited']);

% number of those elements
n_total = length(row);

%% sample event locations
n_sample = 100;

% set seed
s = rng;
p = randperm(n_total);
% location idx
sample_row_idx = row(p(1:n_sample));
sample_col_idx = col(p(1:n_sample));

% event locations - note min_lat -> index 1
sample_lat = min_lat + resol_decimal*(sample_row_idx - 1);
sample_lon = min_lon + resol_decimal*(sample_col_idx - 1);

% check pattern
linear_idx = sub2ind(size(pt_count_mat), sample_row_idx, sample_col_idx);
ind_mat = zeros(size(pt_count_mat));
ind_mat(linear_idx) = 1;
figure;
spy(ind_mat);
set(gca,'YDir','normal');

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

    n_entry(kk) = sum(m_idx);
    
    report(kk).m_idx = m_idx;
    report(kk).locb = locb(m_idx); % the member index wrt event_loc
    report(kk).pt_loc = c_loc(m_idx,:);
    pt_date = data_store(kk).date;
    report(kk).pt_date = pt_date(m_idx);
    pt_time = data_store(kk).time;
    report(kk).pt_time = pt_time(m_idx);
   
%     %% find unique values -- wrong - should be to compress consecutive same values
%     [report(kk).pt_uni_loc,ia,ic] = unique(report(kk).pt_loc, 'rows');
%     % [C,ia,ic]= unique(A)
%     % [~, uni_locb] = ismember(report(kk).pt_uni_loc, event_loc, 'rows');
%     report(kk).uni_locb = report(kk).locb(ia);
%     n_uni_loc(kk) = size(report(kk).pt_uni_loc, 1);
%     report(kk).n_uni_loc = n_uni_loc(kk);
%     report(kk).uni_date = report(kk).pt_date(ia);
%     report(kk).uni_time = report(kk).pt_time(ia);
    
    
    %% trajectory compression
    [report(kk).pt_cpr_loc, idx] = sequence_compress(report(kk).pt_loc);
    report(kk).cpr_locb = report(kk).locb(idx);
    n_cpr_loc(kk) = size(report(kk).pt_cpr_loc, 1);
    report(kk).n_cpr_loc = n_cpr_loc(kk);
    report(kk).cpr_date = report(kk).pt_date(idx);
    report(kk).cpr_time = report(kk).pt_time(idx);
    
    
end

%% visualize how a pt traveled naturally and visited those event locations
scale = 0;
for kk = 1:4 % n_pt
    tmp1 = report(kk).pt_cpr_loc(:,1);
    tmp2 = report(kk).pt_cpr_loc(:,2);
    figure;
    quiver(tmp1, tmp2, [diff(tmp1); 0], [diff(tmp2); 0], scale);
    hold on
    for aa = 1:length(tmp1)
        text(tmp1(aa), tmp2(aa), num2str(report(kk).cpr_locb(aa)));
        text(tmp1(aa), tmp2(aa)+2*resol_decimal, [report(kk).cpr_date(aa) report(kk).cpr_time(aa)]);
    end
    title(['Participant id' num2str(kk)]);
end

