function [tnr, tpr, acc] = binary_acc(conf_mat)

% % calculate true positive rate, true negative rate and accuracy
% neg 1st row and 1st col
% pos 2nd row and 2nd col

n_data = sum(sum(conf_mat));

tnr = conf_mat(1,1)/(conf_mat(1,1) + conf_mat(1,2));
tpr = conf_mat(2,2)/(conf_mat(2,2) + conf_mat(2,1));
acc = (conf_mat(1,1)+conf_mat(2,2))/n_data;
