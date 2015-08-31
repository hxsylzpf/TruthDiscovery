function [class_true, class_est, F1, WAF]= pt_classify_err(d_point, a, a_est, n_pt_class, state)

    class_true = oneD_classify(a, d_point);
    class_est = oneD_classify(a_est, d_point);
    
    [~, confmtx] = confmat(n_pt_class, length(class_true), class_true, state, class_est);
    % note: you cannot use binary here
    [~, ~, F1, WAF, ~] = multiple_f_measure(confmtx);
    