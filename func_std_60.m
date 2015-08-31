function std_comp = func_std_60(data)

% compute the std based on 60% of data

idx = ~isnan(data);
% remove nan
data = data(idx);
% sort data
data = sort(data);

n = length(data);

data_filtered = data(ceil(0.2*n):floor(0.8*n));

std_comp = std(data_filtered);
