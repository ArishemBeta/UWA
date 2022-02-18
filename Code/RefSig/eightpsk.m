function data_sym = eightpsk(data_bit)
%FUNCTION EIGHTPSK
%   data_sym = eightpsk(data_bit)
%

[m, n] = size(data_bit);

map = [3,4,2,1,6,5,7,0];

temp = mod(n, 3);
if temp ~= 0
   for i = 1:3-temp
      n = n+1;
      data_bit(:,n) = 0;
   end
end

for k = 1:n/3
   sym_value(:,k) = 1+data_bit(:, (k-1)*3+3)*4 + data_bit(:, (k-1)*3+2)*2 + data_bit(:,(k-1)*3+1);
end

   sym_value = map(sym_value);
   data_sym = exp(j*pi/4*sym_value);

