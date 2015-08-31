function [data_mat, z] = func_data_mat_filter(data_mat, z, min_count_per_fact, min_count_per_pt)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% events with no label at all (from any user) should not be considered in evaluation
    % label_count_per_fact = sum(~isnan(data_mat));
    
    label_count_per_fact = zeros(1, size(data_mat,2));
    std_per_fact = zeros(1, size(data_mat,2));
    
    for j = 1:size(data_mat,2)
        cur_X = data_mat(:,j);
        label_count_per_fact(j) = sum(~isnan(cur_X));
        % remove corresponding rows immediately
        if label_count_per_fact(j) < min_count_per_fact
           std_per_fact(j) = 0;
        else
           % compute the std of the claims
           std_per_fact(j) = func_std_80(cur_X);
        end
    end
    
    ind_fact_wo_label = (label_count_per_fact < min_count_per_fact) | (std_per_fact < 1);
    
%     % integer index - retained index
%     ind_fact_retain_integer = find(label_count_per_fact >= min_count_per_fact);
    
%     n_fact_wo_label = sum(ind_fact_wo_label);
%     n_s = n_s - n_fact_wo_label;
    z = z(~ind_fact_wo_label);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% update the data matrix
    data_mat = data_mat(:, ~ind_fact_wo_label);

    %% pts with no label at all (for all events) should not be considered as well
    label_count_per_pt = sum(~isnan(data_mat), 2);
    ind_pt_wo_label = (label_count_per_pt < min_count_per_pt); % 0
    
%     % integer index - retained index
%     ind_pt_retain_integer = find(label_count_per_pt >= min_count_per_pt);
    
%     n_pt_wo_label = sum(ind_pt_wo_label);
%     n_pt = n_pt - n_pt_wo_label;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% update the data matrix
    data_mat = data_mat(~ind_pt_wo_label,:);
    