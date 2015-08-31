clear;

file_name = 'data_line';
fileID = fopen([file_name '.csv']);

format_spec = repmat('%q', 1, 8);
C = textscan(fileID,format_spec,'delimiter',',');
fclose(fileID);

worker_id = C{1}(2:end);
line_id = C{4}(2:end);

uni_line_id = unique(line_id);
n_line = length(uni_line_id);

worker_input = C{5}(2:end);
label_num = C{7}(2:end);
label_cat = C{8}(2:end);

% unique labels
uni_state = {'short', 'medium', 'long', 'very long'};
n_state = length(uni_state);

% unique words
uni_word = uni_state;
n_word = length(uni_word);

% unique workers
[uni_worker_id, ia, ic] = unique(worker_id);
n_worker = length(uni_worker_id);

%% filter 1
% filter by num of responses per worker
n_res_per_worker = zeros(n_worker, 1);
for i = 1:n_worker
    n_res_per_worker(i) = sum(strcmp(worker_id, uni_worker_id{i}));
end

idx_retain_worker = (n_res_per_worker >= 5);
retain_uni_worker_id = uni_worker_id(idx_retain_worker);
n_retain_worker = sum(idx_retain_worker)

n_filtered_worker1 = n_worker - n_retain_worker

%% clean all the other data
idx_retain = ismember(worker_id, retain_uni_worker_id);

worker_id = worker_id(idx_retain);
line_id = line_id(idx_retain);

worker_input = worker_input(idx_retain);
label_num = label_num(idx_retain);
label_cat = label_cat(idx_retain);

%% filter 2
% filter by diversity in responses
% number of times a worker uses a word
mat_worker_word = zeros(n_retain_worker, n_word);

for i = 1:n_retain_worker
    for j = 1:n_word
        cur_idx = strcmp(worker_id, retain_uni_worker_id{i}) & strcmp(worker_input, uni_word{j});
        mat_worker_word(i,j) = sum(cur_idx);
    end
end

% filter out workers that use a single kind of word
count_worker_word = sum(mat_worker_word > 0, 2)

% number of times a worker users a num
mat_worker_num = zeros(n_retain_worker, 1);

for i = 1:n_retain_worker
    % find current worker
    cur_idx = strcmp(worker_id, retain_uni_worker_id{i});
    % find current worker's input
    cur_input = worker_input(cur_idx);
    % a num is used if it is not the uni_word set
    cur_idx_num = ~ismember(cur_input, uni_word)
    if sum(cur_idx_num) > 0
        mat_worker_num(i) = 1;
    end
end

% filter out workers that use a single kind of word & no number
filter_uni_worker_id = (count_worker_word == 1) & (mat_worker_num == 0);

%% clean all the other data
idx_retain_worker = ~filter_uni_worker_id;

retain_uni_worker_id = retain_uni_worker_id(idx_retain_worker);

n_retain_worker = sum(idx_retain_worker)

idx_retain = ismember(worker_id, retain_uni_worker_id);

worker_id = worker_id(idx_retain);
line_id = line_id(idx_retain);

worker_input = worker_input(idx_retain);
label_num = label_num(idx_retain);
label_cat = label_cat(idx_retain);

%% num of input
n_input = length(worker_input);
idx_input_word = zeros(n_input, 1);

for i = 1:n_input
    if isempty(str2num(worker_input{i}))
        idx_input_word(i) = true;
    end
end

idx_input_word = logical(idx_input_word);
idx_input_num = ~idx_input_word;

% unique words
% % if you allow words to be arbitrary
% [uni_word, ia_word, ic_word] = unique(worker_input(idx_input_word));
% % if you allow words to be the same as state labels

%% statistics
% how many people use numbers
worker_use_num = unique(worker_id(idx_input_num))
n_worker_use_num = length(worker_use_num)

ratio_worker_use_num = n_worker_use_num/n_retain_worker

% state-dependent number use
% short - more people use numbers
for i = 1:n_state
    cur_state = uni_state{i};
    cur_idx = strcmp(label_cat, cur_state);
    cur_idx_use_num = (cur_idx & idx_input_num);
    cur_user = worker_id(cur_idx_use_num);
    worker_user_num_c_state{i} = unique(cur_user)
    num_c_state{i} = worker_input(cur_idx_use_num)
end

%% user-state-word matrix
count_user_state_word = zeros(n_state, n_word, n_retain_worker);

for i = 1:n_retain_worker
    for j = 1:n_state
        for k = 1:n_word
            cur_idx = strcmp(worker_id, retain_uni_worker_id{i}) & strcmp(label_cat, uni_state{j}) & strcmp(worker_input, uni_word{k});
            cur_count = sum(cur_idx);
            % 1 is a pseudo count
            count_user_state_word(j,k,i) = 0.1 + cur_count;
        end
    end
end

prob_user_state_word = prob_mat_nlz(count_user_state_word, 'row')

%% estimate state log prob
log_prob_line_state = zeros(n_line, n_state);

for i = 1:n_line
    cur_line = uni_line_id{i};
    cur_idx = strcmp(line_id, cur_line);
    
    cur_worker_set = worker_id(cur_idx);
    cur_input_set = worker_input(cur_idx);
    
    n_cur_input = length(cur_worker_set);
    

    for k = 1:n_cur_input
        cur_idx_worker = find(strcmp(cur_worker_set{k}, retain_uni_worker_id))
        cur_idx_word = find(strcmp(cur_input_set{k}, uni_word))
        
        % if not empty
        if ~isempty(cur_idx_word)
           
           for j = 1:n_state
             log_prob_line_state(i,j) = log_prob_line_state(i,j) + log(prob_user_state_word(j,cur_idx_word,cur_idx_worker));
           end
        end
    end
end

%% find max
[max_val, est_state] = max(log_prob_line_state, [], 2)

