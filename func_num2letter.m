function b = func_num2letter(a)
% a is a vec
% length of 24, as required; no J and O
letter = 'ARNDCQEGHILKMFPSTWYVBZX*'; % for 'aa'
% letter = 'ABCDEFGHIKLMNPQRSTUVWXYZ';

b = letter(a);
