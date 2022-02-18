clc;
clear all;
close all;

K=1024;
Nbps=2; %QPSK
Nt=1;
Nr=8;
T214 = poly2trellis(4,[17 13]); %For convolutional encoder
tblen=12;               %For convolutional decoder
rate= 1/2;              %convolutional coding rate
sc_idx= [1:4:K];        %indices of subcarriers for channel estimation
Nbit= K*Nbps*rate;      %number of information bits in one OFDM symbol
SNRdB= 20;              %signal-to-noise ratio in dB
SNR= 10^(SNRdB/10);     %SNR in linear scale
Kg= 240;                %gap between OFDM symbol or between OFDM symbol and pilot block
L= 80;                  %length of channel
Nfrm=1;                % frame
Nblock=1;                  % block per frame
M=4;

%------------generate Tx bits--------------
P_bit=randi([0 1],1,K*Nblock*Nfrm);
% % load P_bit.mat;
bit_seq=reshape(P_bit,K,length(P_bit)/K);
%------------channel coding--------------
for i=1:length(P_bit)/K
code_bitt(i,:)=convenc(bit_seq(:,i),T214);
code_bit(i,:)=randintrlv(code_bitt(i,:),0);
%code_bit=code_bitt;
end
%------------QPSK--------------
data_sym=qpsk(code_bit).';
idx=find(imag(data_sym)==1);
% data_sym=[1;i;-i;1;-1;-i;-1;i];
% for i=1:127
%     data_sym=[data_sym;1;i;-i;1;-1;-i;-1;i];
% end
% scatterplot(data_sym(1:4:end));

%------------IFFT--------------
% ifft_data=ifft(sym_seq); 
ifft_data=ifft(data_sym);

%------------ZP--------------
Tx_block=[ifft_data;zeros(Kg,Nfrm*Nblock)];%zeros(Kg,Nfrm*Nblock);

%------------P/S--------------
Tx_data=(reshape(Tx_block,[],1).');
% plot(abs(Tx_data));

d=0.001;
% Rx_data=interp1((0:length(Tx_data(:))-1),Tx_data(:),(0:length(Tx_data(:))-1)*(1+d),'spline');
Rx_data=resample(Tx_data,1001,1000);
% Rx_data=Tx_data;
Noffset=(1-1)*(K+Kg);%-floor(doppler_scale(max(1,nblk-1))*1264*(nblk-1))
Rx_block=Rx_data(Noffset+1:Noffset+(K+Kg));
% Rx_block=interp1((0:length(Rx_block)-1),Rx_block,(0:length(Rx_block)-1)/(1-d),'spline');
Rx_block=resample(Rx_block,1000,1001);
if(1)
    block_symbol=data_sym.';
    y=Rx_block(:,1: K);
    ola=Rx_block(:,1+K: K+L-1);
    y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
    h= FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, 1, 1, K, L, sc_idx, SNR);
    % plot(abs(h(1,:)));
    LLa_cod= log(0.5)*ones(2*Nt,Nbps*K);
    LLR_info= zeros(Nt,Nbit);
    [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_qpsk_fn_28may11(y,LLa_cod,h,K,1,L,Nt,1,Nbps,SNRdB);
    % scatterplot(S_Est(513:4:640));
    cost1=0;cost2=0;cost3=0;cost4=0;cost5=0;
    symbol_err=S_Est-data_sym.';
    for k=1:256
        cost1=cost1+(symbol_err(k)*conj(symbol_err(k)));
    end
    for k=257:512
        cost2=cost2+(symbol_err(k)*conj(symbol_err(k)));
    end
    for k=513:768
        cost3=cost3+(symbol_err(k)*conj(symbol_err(k)));
    end
    for k=769:1024
        cost4=cost4+(symbol_err(k)*conj(symbol_err(k)));
    end
    for k=385:640
        cost5=cost5+(symbol_err(k)*conj(symbol_err(k)));
    end
    for k=1:K
        cost(k)=(symbol_err(k)*conj(symbol_err(k)));
    end
    plot(cost);
    scatterplot(S_Est(:));
else
    Rx_sym=fft(Rx_block(1:K));
    cost1=0;cost2=0;cost3=0;cost4=0;cost5=0;
    symbol_err=Rx_sym-data_sym.';
%     for k=1:256
%         cost1=cost1+(symbol_err(k)*conj(symbol_err(k)));
%     end
%     for k=257:512
%         cost2=cost2+(symbol_err(k)*conj(symbol_err(k)));
%     end
%     for k=513:768
%         cost3=cost3+(symbol_err(k)*conj(symbol_err(k)));
%     end
%     for k=769:1024
%         cost4=cost4+(symbol_err(k)*conj(symbol_err(k)));
%     end
%     for k=385:640
%         cost5=cost5+(symbol_err(k)*conj(symbol_err(k)));
%     end
    for k=1:K
        cost(k)=(symbol_err(k)*conj(symbol_err(k)));
    end
    plot(cost);
    % for i=1:4:length(Rx_sym)
    %     scatterplot(Rx_sym(i));
    %     hold on;
    % end
    scatterplot(Rx_sym(idx));
    % scatterplot(Rx_sym(129:4:256));
end

% Fs=64;            %采样频率39062.5
% Ts=1/Fs;               %采样间隔
% Fd=1;                  %Doppler shift
% tau=[0.0015,0.003,0.005];          %多径延时向量，s
% pdb=[-10,-18,-22];             %多径信道增益向量，dB
% chan = rayleighchan(Ts, Fd, tau, pdb);
% chan.StoreHistory=1;
% Rxa=filter(chan,a);
% % plot(chan);
% chane=deconv(Rxa,a);

% plot(abs(Rxa));
% hold on;
% Rxa1=interp1([0:length(a)-1],a,[0:length(a)-1]/(1+1),'spline');
% plot(abs(Rxa1));
% for n=1:length(a)
%     Rxa1=a*(exp(sqrt(-1)*2*pi*200*(n-1)/64));
% end
% RxA=fft(Rxa);
% figure();
% plot(abs(RxA));%(1:1/2*length(RxA))

% SNR=[-10,-5,0,5,10,15,20];
% BER=[0.45,0.12,0.002,0,0,0,0];
% figure('color',[1 1 1]);
% plot(SNR,BER,'-*','LineWidth',1.5,'MarkerSize',10);
% ylabel("BER","FontName","Times New Roman","FontSize",20);
% xlabel("SNR(dB)","FontName","Times New Roman","FontSize",20);
% 
% SNR=[-10,-5,0,5,10,15,20];
% RMSE=[6.5635e-4,4.6935e-4,1.5308e-4,1.5264e-4,1.4730e-4,1.4362e-4,1.227e-4];
% figure('color',[1 1 1]);
% plot(SNR,RMSE,'-*','LineWidth',1.5,'MarkerSize',10);
% ylabel("RMSE","FontName","Times New Roman","FontSize",20);
% xlabel("SNR(dB)","FontName","Times New Roman","FontSize",20);

aa=interp1((0:length(seq(:,1))-1),seq(:,1),(0:length(seq(:,1))-1)/(1+0.005),'spline');
bb=xcorr(aa,seq(:,1));plot(bb);