clear;
close all

%% read data
file_name = 'f710250_car_200_copy';
fileID = fopen([file_name '.csv']);

C = textscan(fileID,'%s%s%s%f%f%s','HeaderLines',1,'delimiter',',');
fclose(fileID);

started = C{1};
created = C{2};
worker_id = C{3};
res = C{4};
diff_level = C{5};
task_id = C{6};

uni_worker_id = unique(worker_id)
n_uni_worker = length(uni_worker_id)

% uni_diff_level = unique(diff_level)
% n_uni_diff_level = length(uni_diff_level)

% uni_diff_level = {'Easy', 'Normal', 'Hard'};

uni_diff_level = [1 2 3];

%% construct the data matrix
% X = zeros(n_uni_worker, n_uni_task);
% 
% for i = 1:n_uni_worker
%     for j = 1:n_uni_task
%         idx = (strcmp(worker_id, uni_worker_id{i})) & (strcmp(task_id, uni_task_id{j}));
%         if sum(idx) > 0
%             X(i,j) = str2num(count{idx});
%         end
%     end
% end

%% filter 1
n_min_res = 40;
[n_res_per_worker, idx_retain_worker, idx_retain_res] = num_res_filter(worker_id, uni_worker_id, n_min_res);

my_figure(1/2.5, 1/4);
plot(sort(n_res_per_worker), 'linewidth', 2);
title('Num of res per worker');

n_retain_worker = sum(idx_retain_worker);
n_filtered_worker = n_uni_worker - n_retain_worker

started = started(idx_retain_res);
created = created(idx_retain_res);
worker_id = worker_id(idx_retain_res);
res = res(idx_retain_res);
diff_level = diff_level(idx_retain_res);
task_id = task_id(idx_retain_res);

uni_worker_id = unique(worker_id);
n_uni_worker = length(uni_worker_id);

uni_task_id = unique(task_id);
n_uni_task = length(uni_task_id);

%% reformatting data mat
X_res = nan(n_uni_worker, n_uni_task);
X_diff = nan(n_uni_worker, n_uni_task);

for i = 1:n_uni_worker
    for j = 1:n_uni_task
        cur_idx = (strcmp(worker_id, uni_worker_id{i})) & (strcmp(task_id, uni_task_id{j}));
        if sum(cur_idx) > 0
            X_res(i,j) = res(cur_idx);
            X_diff(i,j) = diff_level(cur_idx);
        end
    end
end

%% statistics 1 - diversity in each user's response
record = struct;
ent_res_whole = zeros(n_uni_worker, 1);
ent_diff_whole = zeros(n_uni_worker, 1);

% for each user, calculate the answer entropy
for i = 1:n_uni_worker
   idx_i = ~isnan(X_res(i,:));
   
   % response
   cur_res = X_res(i,idx_i);
   cur_uni_res = unique(cur_res);
   record(i).uni_res = cur_uni_res;
   record(i).num_res = mycount_unique(cur_res, cur_uni_res);
   % entropy of response
   record(i).prob_res = record(i).num_res/sum(record(i).num_res);
   record(i).ent_res = myentropy_pmf(record(i).prob_res);
   ent_res_whole(i) = record(i).ent_res;
   
   % diff level
   cur_diff = X_diff(i,idx_i);
   record(i).num_diff = mycount_unique(cur_diff, uni_diff_level);
   % entropy of diff level
   record(i).prob_diff = record(i).num_diff/sum(record(i).num_diff);
   record(i).ent_diff = myentropy_pmf(record(i).prob_diff);
   ent_diff_whole(i) = record(i).ent_diff;
   
end

my_figure(1/2.5, 1/4);
plot(sort(ent_res_whole), 'linewidth', 2);
title('Ent res whole');

my_figure(1/2.5, 1/4);
plot(sort(ent_diff_whole), 'linewidth', 2);
title('Ent diff whole');

%% statistics 2 - per task's median value and std
est_z = zeros(n_uni_task, 1);
std_est_z = zeros(n_uni_task, 1);
diff_z = zeros(n_uni_task, 1);
std_diff_z = zeros(n_uni_task, 1);

for j = 1:n_uni_task
    idx_j = ~isnan(X_res(:,j));
    cur_res = X_res(idx_j,j);
    est_z(j) = median(cur_res);
    std_est_z(j) = sqrt(mean((cur_res - est_z(j)).^2));
    
    cur_diff = X_diff(idx_j,j);
    diff_z(j) = mean(cur_diff);
    std_diff_z(j) = sqrt(mean((cur_diff - diff_z(j)).^2));
    
end

my_figure(1/2.5, 1/4);
plot(sort(std_est_z), 'linewidth', 2);
title('Std est z');

my_figure(1/2.5, 1/4);
plot(sort(std_diff_z), 'linewidth', 2);
title('Std diff z');

% correlation btw std_est_z and est_z
% the larger the count, the larger the std since the task is more difficult
std_est_z_nlz = my_min_max_nlz(std_est_z);
est_z_nlz = my_min_max_nlz(est_z);

my_figure(1/2.5, 1/4);
plot(1:n_uni_task, std_est_z_nlz, 1:n_uni_task, est_z_nlz, 'linewidth', 2);
legend('Std est z', 'Est z');
title('correlation btw std_est_z and Est_z');

% correlation btw std_est_z and diff_z
% the larger the perceived diff level, the larger the std
diff_z_nlz = my_min_max_nlz(diff_z);

my_figure(1/2.5, 1/4);
plot(1:n_uni_task, std_est_z_nlz, 1:n_uni_task, diff_z_nlz, 'linewidth', 2);
legend('Std est z', 'Diff z');
title('correlation btw std_est_z and diff_z');

%% statistics 3 - per worker's deviation from median values
std_worker = zeros(n_uni_worker, 1);

for i = 1:n_uni_worker
   idx_i = ~isnan(X_res(i,:));
   
   % response
   cur_res = X_res(i,idx_i);
   
   % diff wrt est_z
   diff_res_est = cur_res - est_z(idx_i).';
   
   % std
   std_worker(i) = sqrt(mean(diff_res_est.^2));
   
end

lambda_worker = 1./std_worker.^2;

my_figure(1/2.5, 1/4);
plot(sort(std_worker), 'linewidth', 2);
title('Std worker');

my_figure(1/2.5, 1/4);
plot(sort(lambda_worker), 'linewidth', 2);
title('Lambda worker');

% new z
z = [ ...
26
32
27
41
36
18
11
10
36
42
53
15
75
17
9
16
12
23
31
21
108
54
34
24
22
19
35
53
17
12
20
136
65
22
38
45
31
47
7
17
10
10
63
15
36
200
56
21
7
37
15
12
60
48
9
12
89
32
10
42
300
7
70
10
17
47
41
6
6
6
15
70
27
25
21
20
72
65
6
31
25
14
35
18
9
54
39
13
22
11
13
37
43
15
39
35
8
27
26
23
6
6
13
30
27
12
123
18
9
16
33
31
12
15
10
200
27
20
12
28
58
23
43
25
10
48
13
10
16
9
71
8
14
55
40
16
104
25
12
45
200
8
42
19
30
17
61
8
10
13
38
15
42
20
136
16
17
19
12
30
25
40
36
105
7
27
15
18
33
12
200
20
9
200
13
35
89
14
77
17
50
24
21
82
102
25
25
33
23
11
90
9
16
23
87
88
15
13
21
87
];


%% remove z>100 - as it goes beyond people's ability
% z = round(est_z);

scale_fac = 1; % 0.1;

idx_retain_s = (z<=100); % (z>= 10 & z <= 60); % 20,50 / 15,100 

z = z(idx_retain_s)*scale_fac;
X_res = X_res(:, idx_retain_s)*scale_fac;
X_diff = X_diff(:, idx_retain_s);
uni_task_id = uni_task_id(idx_retain_s);

%% calculate statistics 3 again - per worker's deviation from true value
std_worker = zeros(n_uni_worker, 1);

for i = 1:n_uni_worker
   idx_i = ~isnan(X_res(i,:));
   
   % response
   cur_res = X_res(i,idx_i);
   
   % diff wrt est_z
   diff_res_est = cur_res - z(idx_i).';
   
   % std
   std_worker(i) = sqrt(mean(diff_res_est.^2));
   
end

lambda_worker = 1./std_worker.^2;

my_figure(1/2.5, 1/4);
plot(sort(std_worker), 'linewidth', 2);
title('Std worker new');

my_figure(1/2.5, 1/4);
plot(sort(lambda_worker), 'linewidth', 2);
title('Lambda worker new');

%% save data
save('data_car_cf_all_20150420.mat', 'X_res', 'X_diff', 'z', 'uni_task_id', 'uni_worker_id', 'scale_fac');
