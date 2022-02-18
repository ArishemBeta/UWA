SampleRate=9765.625;            %采样频率39062.5
MaximumDopplerShift=10;                  %Doppler shift
NumTransmitAntennas=Nt;
NumReceiveAntennas=Nr;
RandomStream='mt19937ar with seed';
PathDelays=[0.0015,0.003];          %多径延时向量，s
AveragePathGains=[-10,-18];             %多径信道增益向量，dB
mimochan = comm.MIMOChannel('SampleRate',SampleRate,'MaximumDopplerShift',MaximumDopplerShift,'PathDelays',PathDelays,...
'AveragePathGains',AveragePathGains,'NumTransmitAntennas',NumTransmitAntennas,'NumReceiveAntennas',NumReceiveAntennas,...
'RandomStream','Global stream','SpatialCorrelationSpecification','None');%mt19937ar with seed
Rx_data=mimochan(Tx_data).';
% Rx_data=Tx_data;,'DopplerSpectrum',doppler('Flat')
Rx_data=awgn(Rx_data,SNRdB,'measured');