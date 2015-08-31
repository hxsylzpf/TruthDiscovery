function merged_diff = func_merge_level(per_diff, style)

% merge 12 as 1 level
if strcmp(style, '12,3')
   merged_diff = per_diff;
   idx1 = (per_diff == 1);
   idx2 = (per_diff == 2);
   idx3 = (per_diff == 3);
   
   idx12 = idx1 | idx2;
   merged_diff(idx12) = 1;
   merged_diff(idx3) = 2;

% merge 23 as 1 level   
elseif strcmp(style, '1,23')
   merged_diff = per_diff;
   idx1 = (per_diff == 1);
   idx2 = (per_diff == 2);
   idx3 = (per_diff == 3);
   
   idx23 = idx2 | idx3;
   merged_diff(idx1) = 1;
   merged_diff(idx23) = 2;
    
end