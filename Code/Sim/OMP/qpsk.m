function data_sym = qpsk(data_bit)
%FUNCTION QPSK
%   data_sym = qpsk(data_bit)
%

if and(data_bit(1) ~= 0, data_bit(1) ~= 1)
   if and(data_bit(1) ~= '0', data_bit(1) ~= '1')
      error('data_bit must be either integer (1, 0), or ascii sym (1, 0)');
   else
      data_bit = data_bit-48;
   end
end

  

[m, n] = size(data_bit);

map = [0, 1, 3, 2];

temp = mod(n, 2);
if temp ~= 0
   n = n+1;
   data_bit(:,n) = 0;
end

sym_value = zeros(m, n/2);
for k = 1:n/2
   sym_value(:,k) = 1+data_bit(:, (k-1)*2+1)*2 + data_bit(:,(k-1)*2+2);
end

sym_value = map(sym_value);
data_sym = exp(j*pi/2*sym_value);

