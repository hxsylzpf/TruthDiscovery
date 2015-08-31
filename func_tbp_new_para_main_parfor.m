function record = func_tbp_new_para_main_parfor(files, n_node, n_pt, n_l, ini_pi, ini_h, ini_lambda, prior_a, prior_b)

% n_node is the number of nodes to process the data in parallel

n_iter = 100;

pi = ini_pi;
h = ini_h;
lambda = ini_lambda;

map_time = 0;
reduce_time = 0;
check_time = 0;

for kkk = 1:n_iter
    %% mapper
    tic
    func_tbp_new_para_mapper_parfor(files, n_l, pi, h, lambda);
    temp_time1 = toc;
    
    map_time = map_time + temp_time1;
    
    %% read the mapper output
    mapper_files = cell(n_node,1);

    for tt = 1:n_node
        mapper_files{tt} = ['mapper_output_split_' num2str(n_node) '_file_' num2str(tt) '.mat'];
    end
    
    %% reducer
    tic
    reducer_record = func_tbp_new_para_reducer_parfor(mapper_files, n_node, n_pt, n_l, prior_a, prior_b);
    temp_time2 = toc;
    
    reduce_time = reduce_time + temp_time2;
    
    tic
    new_pi = reducer_record.pi_est;
    new_h = reducer_record.h_est;
    new_lambda = reducer_record.lambda_est;
    z_est = reducer_record.z_est;
    
    %% check convergence
    if norm(new_h(:) - h(:)) + norm(new_lambda(:) - lambda(:)) < 1e-3
       break;
    else
       pi = new_pi;
       h = new_h;
       lambda = new_lambda;
    end
        
    temp_time3 = toc;
    
    check_time = check_time + temp_time3;
end

cpu_time = map_time + reduce_time + check_time;

% record.pi_est = pi;
% record.h_est = h;
% record.lambda_est = lambda;

record.map_time = map_time;
record.reduce_time = reduce_time;
record.cpu_time = cpu_time;
record.z_est = z_est;
