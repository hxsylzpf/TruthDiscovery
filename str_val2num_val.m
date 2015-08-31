function num_val = str_val2num_val(str_val, str_cat, num_cat)
% map string values e.g., 'A', 'B' to specified num values e.g., '1', '2'
% str_cat - string category
% num_cat - num category

% num of label category
n_cat = length(num_cat);
% 
% if n_cat~=length(uni_str_val)
%     error('Error - num of unique labels do not equal num of categories');
% end   

num_val = zeros(size(str_val));

for i = 1:n_cat
    cur_idx = strcmp(str_val, str_cat{i});
    num_val(cur_idx) = num_cat(i);
end
