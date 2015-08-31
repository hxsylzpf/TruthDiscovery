function std_comp = func_std_80(data)

% compute the std based on 80% of data

idx = ~isnan(data);
% remove nan
data = data(idx);
% sort data
data = sort(data);

n = length(data);

data_filtered = data(ceil(0.1*n):floor(0.9*n));

std_comp = std(data_filtered);
