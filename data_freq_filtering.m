clear;

% filter out pts based on freq of data

%% load data set
dataset_name = 'data_count_cf_15_100_20141006'; % ''; % 'data_count_cf_new_15_100_20150128'; % 'data_car_cf_all_20150420';
a = load([dataset_name '.mat']);
X = a.X_res;
per_diff = a.X_diff;

n_pt = size(X,1);

n_l = 3;

n_claim_per_level = zeros(n_pt,n_l);

% %% remove 0 in data
% idx0 = (X == 0);
% X(idx0) = nan;
% per_diff(idx0) = nan;

%% setting
% two different settings of min_claim_per_level
min_claim_per_level = 5;
min_claim_per_level2 = 10;

% copy ori
X_filtered = X;
diff_filtered = per_diff;
n_claim_per_level_filtered = zeros(n_pt,n_l);

X_filtered2 = X;
diff_filtered2 = per_diff;
n_claim_per_level_filtered2 = zeros(n_pt,n_l);

%% compute the num of claims per diff level
for i = 1:n_pt
    for k = 1:n_l
        idx = (per_diff(i,:) == k);
        n_claim_per_level(i,k) = sum(idx);
        
        % if smaller than a threshold, set X and diff to nan
        if n_claim_per_level(i,k) < min_claim_per_level
           X_filtered(i,idx) = nan;
           diff_filtered(i,idx) = nan;
        end
        
        if n_claim_per_level(i,k) < min_claim_per_level2
           X_filtered2(i,idx) = nan;
           diff_filtered2(i,idx) = nan;
        end
        
    end
end

%% recompute the num of claims per diff level - check correctness
for i = 1:n_pt
    for k = 1:n_l
        idx = (diff_filtered(i,:) == k);
        n_claim_per_level_filtered(i,k) = sum(idx);
    end
end

%% remove pts whose diff level is not compelte
idx_pt = (sum(n_claim_per_level_filtered == 0, 2) > 0)

X_filtered(idx_pt,:) = [];
diff_filtered(idx_pt,:) = [];

% % if car
% X_filtered(end,:) = [];
% diff_filtered(end,:) = [];

save(dataset_name, 'X_filtered', 'diff_filtered', '-append'); % 'X_filtered2', 'diff_filtered2', 
