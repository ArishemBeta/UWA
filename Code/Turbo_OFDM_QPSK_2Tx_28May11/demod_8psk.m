function data = demod_8psk(sym)
%FUNCTION DEMOD_8PSK
%   data = demod_8psk(sym)
%   demodulate 8psk sym, and return the binary data
%

[rows, columns] = size(sym);

%map = [7,3,2,0,1,5,4,6];
map = [7, 6, 2, 0, 4, 5, 1, 3];

phase_sym = angle(sym); %this angle is between [-pi, pi]
%to make the angle between [0, 2*pi]
temp = 2*pi*(phase_sym < 0);
phase_sym = phase_sym+temp;

Oct_matrix = zeros(rows, columns);
for m = 2:8
   Oct_matrix = Oct_matrix+(m-1)*and(phase_sym>(m-1)*pi/4-pi/8, phase_sym<=(m-1)*pi/4+pi/8);
end

Oct_matrix = Oct_matrix+1;
Oct_matrix = map(Oct_matrix);

for m = 1:rows
   temp = dec2bin(Oct_matrix(m,:),3);
   bin_matrix(m,1:3) = temp(1,:);
   for n = 2:columns
      bin_matrix(m,1:n*3) = [bin_matrix(m,1:(n-1)*3), temp(n,1:3)];
   end
end

data = bin_matrix-48;

