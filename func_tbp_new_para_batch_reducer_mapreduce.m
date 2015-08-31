function func_tbp_new_para_batch_reducer_mapreduce(intermKey, intermValIter, outKVStore, n_l)

% pi
if strcmp(intermKey, 's0')

    sum_val = zeros(1, n_l+1);
    while hasnext(intermValIter)
        sum_val = sum_val + getnext(intermValIter);
    end
    
    pi_val = (sum_val(1:n_l) + 1)/(sum_val(end) + n_l);
    add(outKVStore, 'pi', pi_val);
    
% h and lam
else
    sum_val = zeros(1, 3*n_l);

    while hasnext(intermValIter)
        sum_val = sum_val + getnext(intermValIter);
    end
        
    s1 = sum_val(1:n_l);
    s2 = sum_val(n_l+1:2*n_l);
    s3 = sum_val(2*n_l+1:3*n_l);
    
    h = s2./s1;
    lambda = (0.5*s1 + 1e-4)./(0.5*s1.*h.^2 - s2.*h + 0.5*s3 + 1e-4);
    
    add(outKVStore, ['h' intermKey], h);
    add(outKVStore, ['lam' intermKey], lambda);
end
