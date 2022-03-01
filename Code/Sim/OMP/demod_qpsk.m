function data = demod_qpsk(sym)
%FUNCTION DEMOD_8PSK
%   data = demod_qpsk(sym)
%   demodulate qpsk sym, and return the binary data
%

[rows, columns] = size(sym);

map = [0,1,3,2];

phase_sym = angle(sym); %this angle is between [-pi, pi]
%to make the angle between [0, 2*pi]
temp = 2*pi*(phase_sym < 0);
phase_sym = phase_sym+temp;

Oct_matrix = zeros(rows, columns);
for m = 2:4
   Oct_matrix = Oct_matrix+(m-1)*and(phase_sym>(m-1)*pi/2-pi/4, phase_sym<=(m-1)*pi/2+pi/4);
end

Oct_matrix = Oct_matrix+1;
Oct_matrix = map(Oct_matrix);

for m = 1:rows
   temp = dec2bin(Oct_matrix(m,:),2)-48;
   bin_matrix(m,1:2) = temp(1,:);
   for n = 2:columns
      bin_matrix(m,1:n*2) = [bin_matrix(m,1:(n-1)*2), temp(n,1:2)];
   end
end

data = bin_matrix;

