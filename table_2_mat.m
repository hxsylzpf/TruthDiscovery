function data_mat = table_2_mat(data_table)

% pt_id, fact_id must be numeric and sequentially numbered

% [pt_id, fact_id, val]
% convert a data table to an equivalent mat format

% uni_pt_id = unique(data_table(:,1));
% uni_fact_id = unique(data_table(:,2));
% 
% n_pt = length(uni_pt_id);
% n_fact = length(uni_fact_id);

max_n_pt = max(data_table(:,1));
max_n_fact = max(data_table(:,2));

data_mat = nan(max_n_pt, max_n_fact);

n_data = size(data_table,1);

for i = 1:n_data
%     new_pt_id = (data_table(i,1) == uni_pt_id);
%     new_fact_id = (data_table(i,2) == uni_fact_id);
%     data_mat(new_pt_id, new_fact_id) = data_table(i,3);

data_mat(data_table(i,1), data_table(i,2)) = data_table(i,3); 
end

