clear;

%% you need to specify the file name
filename = 'f384204.csv';
sheet = 1;
%% you need to change here (range to read) as well
xlRange = 'I2:R323';

[num,txt,raw] = xlsread(filename, sheet, xlRange);

% get data
pt_ID = num(:,2);
pt_country = raw(:,4);
pt_region = raw(:,5);
pt_label = raw(:,8);
event_ID = raw(:,10);

% convert str category labels to num labels - cate1 -> 1, cate2 -> 0
num_pt_label = str_val2num_val(pt_label, [1 0]);
[n_pt, n_event, uni_pt_ID, uni_event_ID, data_mat, loc_cover_mat] = data_loc_cover_mat_gen(pt_ID, event_ID, num_pt_label);

% loc_cover_mat(8,:) = loc_cover_mat(8,:) | loc_cover_mat(14,:);
% data_mat(8,:) = data_mat(8,:) | data_mat(14,:);
% loc_cover_mat([4 14],:) = [];
% data_mat([4 14],:) = [];
% n_pt = n_pt - 2;

% event label needs to be manually provided
% % task 1: bike rack
% event_label = [ ...
%     1 0 0 0 1 1 1 1 1 0 ...
%     1 0 1 1 0 0 0 1 0 0 ...
%     1 1 0 1 0 1 0 0 1 1];

% task 2: cherry label
event_label = [ ...
    1 1 0 0 1 1 0 0 1 0 ...
    0 1 1 0 1 0 1 1 1 0 ...
    1 1 0 0 0 1 0 1 1 0];

save('CF_task1_2.mat', 'uni_pt_ID', 'uni_event_ID', 'n_pt', 'n_event', 'data_mat', 'loc_cover_mat', 'event_label');
