function class = oneD_classify(input, d_point)
% input - one D vector of nums
% d_point - the points to differentiate classes - btw 0 and 1
n_class = length(d_point) + 1;
n_inst = length(input);
class = zeros(n_inst, 1);
range = [0 d_point 1];

for i = 1:n_class
    idx = (input >= range(i)) & (input < range(i+1));
    class(idx) = i;
    
    % 1 is not considered in the above calculation
    idx1 = (input == 1);
    class(idx1) = n_class;
end
