function data_sym = bpsk(data_bit)
%FUNCTION BPSK
%   data_sym = bpsk(data_bit)
%

[m,n] = size(data_bit);

data_sym = ones(m, n);
data_sym(data_bit == 0) = -1; 

