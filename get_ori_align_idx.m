function ori_align_idx = get_ori_align_idx(align_reveal_seq, align_sign)
% find the original align index for the two aligned sequences

        % relabel and ignore gaps
        len_align = length(align_reveal_seq);
        ori_align_idx = [];
        pointer = 0;
        for i = 1:len_align
            cur_align_reveal = align_reveal_seq(i);
            cur_sign = align_sign(i);
            
            % only when it is not -, move the index for one
            if cur_align_reveal ~= '-'
                pointer = pointer + 1;
            end
            
            if cur_sign == '|' && pointer ~= 0
                ori_align_idx = [ori_align_idx pointer];
            end
        end