clc;
clear all;
close all;

K=1024;
B=9765.625;             %带宽
Ns=8;                   %过采样倍数
Kg= 240;                %保护间隔
Nblock=2;               %OFDM块数
Nfrm=1;
%----随机生成QPSK符号----
data_sym_t=zeros(K,Nblock);
temp=zeros(K,Nblock);
for i=1:Nblock
    for j=1:K
        temp(j,i)=randi([1 4]);
        if(temp(j,i)==1) data_sym_t(j,i)=1+0*sqrt(-1);
        elseif(temp(j,i)==2) data_sym_t(j,i)=0+1*sqrt(-1);
        elseif(temp(j,i)==3) data_sym_t(j,i)=-1+0*sqrt(-1);
        else data_sym_t(j,i)=0+-1*sqrt(-1);
        end
    end
end
% scatterplot(data_sym_t(:,1));
% for i=1:16
%     data_sym_t(K/2+i,:)=0;
%     data_sym_t(K/2+1-i,:)=0;
%     data_sym_t(K/4+i,:)=0;
% end

%------------频域插零----------
% data_sym=data_sym_t;
data_sym=zeros(K*Ns,Nblock);
for i=1:Nblock                              %在中间插零，使ifft后实现Ns倍时域过采样
    for j=1:K/2
        data_sym(j,i)=data_sym_t(j,i);
    end
    for j=K/2+1:K
        data_sym(j+(Ns-1)*K,i)=data_sym_t(j,i);
    end
end

%------------ifft----------
ifft_data=ifft(data_sym);

%------------插入保护间隔----------
Tx_block=[ifft_data;zeros(Kg*Ns,Nfrm*Nblock)];

Tx_data_t=reshape(Tx_block,[],1).';

%-----------调制到通带----------
for i=1:length(Tx_data_t)
    Tx_data(i)=real(Tx_data_t(i)*(exp(sqrt(-1)*2*pi*13000*(i-1)/(B*Ns))));%采样间隔为1/(B*Ns)
end
% fft(Tx_data(1:Ns:K*Ns));

% plot(Tx_data);
%------------信道----------
% Npath=1;
% UWAchannel=UWAchannel_generation(Npath,Nblock*(K+Kg)/B,1/(Ns*B),80,0.1,4,1);
% Rx_data=zeros(1,Nblock*(K+Kg)*Ns+L*Ns);
% for t=1:height(UWAchannel)
%     for p=1:Npath
%         if(t-UWAchannel(t,p*2)*Ns<0)
%             Rx_data(t)=0;
%         else
%             Rx_data(t)=Rx_data(t)+UWAchannel(t,p*2-1)*Tx_data(round(t-UWAchannel(t,p*2)*Ns));
%         end
%     end
% end
% Rx_data=awgn(Tx_data,SNRdB,'measured');
Rx_data=Tx_data;
% figure();
% plot(Rx_data);

%------------转基带-----------
% Rx_data=Rx_data+sqrt(-1)*Rx_data_i;
for i=1:length(Rx_data)
    Rx_data(i)=Rx_data(i)*(exp(sqrt(-1)*2*pi*(-13000)*(i-1)/(B*Ns)));%
end
% Rx_data_i=hilbert(Rx_data);
Rx_data=lowpass(Rx_data,B*0.8,B*Ns);

for i=1:1%Nblock
    Rx_block=Rx_data(1:Ns:K*Ns);
    Rx_sym=fft(Rx_block);
    bit=demod_qpsk(data_sym_t(:,i).').';
    Rx_bit=demod_qpsk(Rx_sym).';
    err=sum(xor(bit,Rx_bit));
    scatterplot(Rx_sym);
end

