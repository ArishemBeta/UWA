SampleRate=9765.625;            %����Ƶ��39062.5
MaximumDopplerShift=10;                  %Doppler shift
NumTransmitAntennas=Nt;
NumReceiveAntennas=Nr;
RandomStream='mt19937ar with seed';
PathDelays=[0.0015,0.003];          %�ྶ��ʱ������s
AveragePathGains=[-10,-18];             %�ྶ�ŵ�����������dB
mimochan = comm.MIMOChannel('SampleRate',SampleRate,'MaximumDopplerShift',MaximumDopplerShift,'PathDelays',PathDelays,...
'AveragePathGains',AveragePathGains,'NumTransmitAntennas',NumTransmitAntennas,'NumReceiveAntennas',NumReceiveAntennas,...
'RandomStream','Global stream','SpatialCorrelationSpecification','None');%mt19937ar with seed
Rx_data=mimochan(Tx_data).';
% Rx_data=Tx_data;,'DopplerSpectrum',doppler('Flat')
Rx_data=awgn(Rx_data,SNRdB,'measured');