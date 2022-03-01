function data_mod = mqam(data, M)
% function data_mod = mqam(data, M, Es)
%    M: 16 or 64
%

Es = 1; % the energy per symbol
% the size of one decision area 2*delta x 2*delta
delta = sqrt(3*Es/(2*(M-1)));

N_data = length(data);


if M == 4
   mapping = [-1, 1]*delta;

   data_inphase = data(1:2:end);
   data_sym_inphase = mapping(data_inphase+1);

   data_quadrature = data(2:2:end);
   data_sym_quadrature = mapping(data_quadrature+1);
elseif M == 16
   mapping = [-3, -1, 3, 1]*delta;

   data_inphase_1 = data(1:4:end);
   data_inphase_2 = data(2:4:end);
   data_sym_inphase = mapping((data_inphase_2 + data_inphase_1*2)+1);
   
   data_quadrature_1 = data(3:4:end);
   data_quadrature_2 = data(4:4:end);
   data_sym_quadrature = mapping((data_quadrature_2 + data_quadrature_1*2)+1);
elseif M == 64
   mapping = [-7, -5, -1, -3, 7, 5, 1, 3]*delta;
   
   data_inphase_1 = data(1:6:end);
   data_inphase_2 = data(2:6:end);
   data_inphase_3 = data(3:6:end);
   data_sym_inphase = mapping((data_inphase_3+data_inphase_2*2+data_inphase_1*4)+1);
   
   data_quadrature_1 = data(4:6:end);
   data_quadrature_2 = data(5:6:end);
   data_quadrature_3 = data(6:6:end);
   data_sym_quadrature = mapping((data_quadrature_3+data_quadrature_2*2+data_quadrature_1*4)+1);
end

data_mod = data_sym_inphase+j*data_sym_quadrature;
  
