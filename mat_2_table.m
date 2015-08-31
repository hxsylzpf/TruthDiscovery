function data_table = mat_2_table(data_mat)

% convert a data matrix to an equivalent table format
[n_pt, n_fact] = size(data_mat);
% n_data = sum(~isnan(data_mat(:)));

% [pt_id, fact_id, val]
data_table = [];
h = 0;

% data table will be sorted in terms of fact; not pt
% that is why fact should be looped outside
for j = 1:n_fact
    for i = 1:n_pt
        % if not nan
        if ~isnan(data_mat(i,j))
            h = h+1;
            data_table(h,:) = [i j data_mat(i,j)];
        end
    end
end

