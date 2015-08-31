function [data_table, z] = func_data_table_filter(data_table, z, min_count_per_fact, min_count_per_pt)

%% calculate stat for facts and pts
% max pt and s
max_n_pt = max(data_table(:,1));
max_n_s = max(data_table(:,2));

%% events with no label at all (from any user) should not be considered in evaluation
% events with 0 std should not be considered in evaluation - there is no need for truth discovery

label_count_per_fact = zeros(1, max_n_s);
std_per_fact = zeros(1, max_n_s);
for j = 1:max_n_s
    idx_j = (data_table(:,2) == j);
    label_count_per_fact(j) = sum(idx_j);
    % remove corresponding rows immediately
    if label_count_per_fact(j) < min_count_per_fact
       data_table(idx_j,:) = [];
       std_per_fact(j) = 0;
    else
       % compute the std of the claims
       std_per_fact(j) = func_std_80(data_table(idx_j,3));
    end
end

%%
ind_fact_wo_label = (label_count_per_fact < min_count_per_fact) | (std_per_fact < 1);
% set corresponding z to nan - not to consider
z(ind_fact_wo_label) = nan;

%% remove the last few Nans, as they will not be considered by other methods
n_s_retained = max(data_table(:,2));
z = z(1:n_s_retained);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pts with no label at all (for all events) should not be considered as well
label_count_per_pt = zeros(max_n_pt, 1);
for i = 1:max_n_pt
    idx_i = (data_table(:,1) == i);
    label_count_per_pt(i) = sum(idx_i);
    % remove corresponding rows immediately
    if label_count_per_pt(i) < min_count_per_pt
       data_table(idx_i,:) = [];
    end
end
