function record = func_tbp_new_para_batch_mapreduce(data_file_name, n_pt, n_l, ini_pi, ini_h, ini_lambda)

n_iter = 50;

%% data read
ds = datastore(data_file_name);

pi = ini_pi;
h = ini_h;
lambda = ini_lambda;

% p = parpool('local', 2);
% mr = mapreducer(p);
% mapreduce(x,x,mr)

tic
for kkk = 1:n_iter
%% map reduce
result = mapreduce(ds, @(data,~,intermKVStore)func_tbp_new_para_batch_mapper_mapreduce(data, [], intermKVStore, n_l, pi, h, lambda), ...
                       @(intermKey, intermValIter, outKVStore)func_tbp_new_para_batch_reducer_mapreduce(intermKey, intermValIter, outKVStore, n_l));

result_tab = readall(result);

result_key = result_tab.Key;
result_val = result_tab.Value;

%% assign vals to vars
% sort_key is used for debug
[sort_key, idx] = sort(result_key);

sort_val = result_val(idx);

sort_val = cell2mat(sort_val);

new_h = sort_val(1:n_pt,:);
new_lambda = sort_val(n_pt+1:2*n_pt,:);
new_pi = sort_val(end,:);
    
    %% check convergence
    if norm(new_h(:) - h(:)) + norm(new_lambda(:) - lambda(:)) < 1e-3
       break;
    else
       pi = new_pi;
       h = new_h;
       lambda = new_lambda;
    end

end
cpu_time = toc;

% p = gcp;
% delete(p);

record.h = h;
record.lambda = lambda;
record.pi = pi;
record.cpu_time = cpu_time;

