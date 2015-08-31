function result = per_loc_tend(n_bin, n_pt, n_event, event_loc_xy, provided_label_mat, loc_cover_mat)

loc_cen_record = zeros(n_pt, 2);
p_r_record = zeros(n_pt, n_bin);
r_dist_record = cell(n_pt, 1);
v_dist_record = cell(n_pt, 1);
r_dist_record_all = [];
v_dist_record_all = [];

personal_loc_tend_mat = zeros(n_pt, n_event);

for i = 1:n_pt
    loc_cen_record(i,:) = median(event_loc_xy(logical(provided_label_mat(i,:)),:), 1);
    
%     % examine revealed locs and loc center
%     figure;
%     hold on
%     plot(event_loc_xy(logical(provided_label_mat(i,:)),1), event_loc_xy(logical(provided_label_mat(i,:)),2), 'bo');
%     plot(loc_cen_record(i,1), loc_cen_record(i,2), 'rd');
%     hold off

    r_dist_record{i} = pdist2(loc_cen_record(i,:), event_loc_xy(logical(provided_label_mat(i,:)),:), 'cityblock');
    v_dist_record{i} = pdist2(loc_cen_record(i,:), event_loc_xy(logical(loc_cover_mat(i,:)),:), 'cityblock');
    
    r_dist_record_all = [r_dist_record_all r_dist_record{i}];
    v_dist_record_all = [v_dist_record_all v_dist_record{i}];
    
%     figure;
%     hist(r_dist_record{i});
%     figure;
%     hist(v_dist_record{i});
    
    % hist
    [nele_r,xcen_r] = hist(r_dist_record{i}, n_bin);
    % prob
    p_r_record(i,:) = (nele_r+0.01)/sum(nele_r+0.01);
    
    % the dist to predict
    dist_to_pre = pdist2(loc_cen_record(i,:), event_loc_xy, 'cityblock');
    % find bin dist
    bin_dist = pdist2(dist_to_pre.', xcen_r.');
    [min_val, bin_idx] = min(bin_dist, [], 2);
    % check the prob table and determine f_ij
    % scaling
    temp = p_r_record(i,:);
    if max(temp)~= min(temp)
        p_r_record(i,:) = 0.1 + 0.8*(temp - min(temp))/(max(temp) - min(temp));
    end
    personal_loc_tend_mat(i,:) = p_r_record(i,bin_idx);   
end

result.personal_loc_tend_mat = personal_loc_tend_mat;
result.loc_cen_record = loc_cen_record;
result.r_dist_record = r_dist_record;
result.v_dist_record = v_dist_record;
result.r_dist_record_all = r_dist_record_all;
result.v_dist_record_all = v_dist_record_all;
result.p_r_record = p_r_record;
