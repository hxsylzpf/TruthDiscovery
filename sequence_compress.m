function [c_seq, idx] = sequence_compress(seq)
% compress a sequence like [1 1 1 2 2 3 3 3]
% as [1 2 3] with indices [1 4 6]
% how many rows
n = size(seq, 1);
if n == 0
    c_seq = [];
    idx = [];   
else
    c_seq = seq(1,:);
    idx = 1;

    for i = 1:n-1
        tmp1 = seq(i,:);
        tmp2 = seq(i+1,:);
        if sum(tmp2 - tmp1) ~=0
            c_seq = [c_seq; tmp2];
            idx = [idx; i+1];
        end
    end
end


    