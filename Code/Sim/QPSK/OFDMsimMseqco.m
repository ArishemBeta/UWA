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
sc_idx= [1:4:K];        %indices of subcarriers for channel estimation [1:2:K/2,K/2+1:8:K]
Nbit= K*Nbps*rate;      %number of information bits in one OFDM symbol
SNRdB= 10;              %signal-to-noise ratio in dB
SNR= 10^(SNRdB/10);     %SNR in linear scale
Kg= 265;                %gap between OFDM symbol or between OFDM symbol and pilot block 127/128 256/265 512/539
L= 80;                  %length of channel
Nfrm=1;                % frame
Nblock=10;                  % block per frame
M=4;

%------------generate Tx bits--------------
P_bit=randi([0 1],1,K*Nblock*Nfrm);
% load P_bitQPSK.mat;
bit_seq=reshape(P_bit,K,length(P_bit)/K);

%------------channel coding--------------
for i=1:length(P_bit)/K
code_bitt(i,:)=convenc(bit_seq(:,i),T214);
code_bit(i,:)=randintrlv(code_bitt(i,:),0);
%code_bit=code_bitt;
end

%------------QPSK--------------
data_sym=qpsk(code_bit).';

%------------S/P--------------
% sym_seq=reshape(data_sym,K,length(data_sym)/K);

%------------IFFT--------------
ifft_data=ifft(data_sym);

%------------ZP--------------
Tx_block=[ifft_data;zeros(Kg,Nfrm*Nblock)];%zeros(Kg,Nfrm*Nblock);
% coef=dec2bin(primpoly(7));
% for i=1:127
%     seq(i)=abs(sin(2*pi*i/256));
% end
% seq=seq.';
seq=mseq([1 0 0 0 1 1 1 0]).';%1 0 0 0 0 1 0 0 0       1 1 0 0 1 0 1
seq=0.06*[seq,seq,seq,seq,seq,seq,seq,seq,seq,seq];
% pb=randi([0 1],1,K*2);
% ps=qpsk(pb).';
% seq=ifft(ps);
% seq=1*[seq,seq,seq,seq,seq,seq,seq,seq,seq,seq];
Tx_block=[seq;zeros(Kg,Nfrm*Nblock);Tx_block];
%------------P/S--------------
Tx_data=reshape(Tx_block,[],1).';
Tx_data=[Tx_data,zeros(1,Kg)];
Tx_data=(Tx_data);
% figure(1);
% plot(abs(Tx_data));

%------------Rayleigh channel--------------
SampleRate=9765.625;            %采样频率39062.5
MaximumDopplerShift=0.1;                  %Doppler shift
NumTransmitAntennas=Nt;
NumReceiveAntennas=Nr;
RandomStream='mt19937ar with seed';
PathDelays=[0.000,0.0008,0.0015];          %多径延时向量,0.008
AveragePathGains=[0,-3,-6];             %多径信道增益向量,-20
mimochan = comm.MIMOChannel('SampleRate',SampleRate,'MaximumDopplerShift',MaximumDopplerShift,'PathDelays',PathDelays,...
'AveragePathGains',AveragePathGains,'NumTransmitAntennas',NumTransmitAntennas,'NumReceiveAntennas',NumReceiveAntennas,...
'RandomStream','mt19937ar with seed','Seed',8001,'SpatialCorrelationSpecification','None');%Global stream
% for n=1:length(Tx_data)
%     Tx_data(n)=(Tx_data(n)*(exp(sqrt(-1)*2*pi*13000*(n-1)/9765.625)));%48828.125
% end

Rx_datat=mimochan(Tx_data.').';
% Rx_datat=Tx_data.';
d=0.005;%Nr=1;
for i=1:Nr
%     Rx_data(i,:)=interp1((0:length(Rx_datat(i,:))-1),Rx_datat(i,:),(0:length(Rx_datat(i,:))-1)*(1+d),'spline');
    Rx_data(i,:)=resample(Rx_datat(i,:),1000,1000+1000*d);
end
Rx_data=awgn(Rx_data,SNRdB,'measured');
c=-5;
% for i=1:Nr
%     for n=1:length(Rx_data(i,:))
%         Rx_data(i,n)=Rx_data(i,n)*(exp(sqrt(-1)*2*pi*(-13000)*(n-1)/9765.625));%48828.125
%     end
% end
% Rx_data=Tx_data;
% plot(abs(Rx_data(1,:)));
%------------S/P--------------
doppler_scale=zeros(1,Nblock);
cfo=zeros(1,Nblock);
shift=zeros(1,Nblock);
% doppler_scale(1)=0.005;
for nblk=1: Nblock %1
    nblk
    Noffset=(nblk-1)*(K+2*Kg+length(seq));
%     Noffset=Noffset-floor(doppler_scale(1)*(K+2*Kg+length(seq)));
    if(nblk>1)
        for i=2:nblk
            Noffset=Noffset-round(doppler_scale(i-1)*(K+2*Kg+length(seq)));%start of the nblk-th block
        end
        Rx_blockt=Rx_data(:,Noffset+1:Noffset+round((K+Kg+Kg+length(seq))*(1-doppler_scale(nblk-1))));
    else
        Rx_blockt=Rx_data(:,Noffset+1:Noffset+(K+Kg+Kg+length(seq)));
    end
    block_symbol=data_sym(:,nblk).';
    pilot_symbol=block_symbol(sc_idx);
    block_bit=bit_seq(:,nblk).';
%     pilot_bit=block_bit(sc_idx);
    
%     plot(abs(Rx_blockt(1,:)));
    [doppler_scale(nblk),shift(nblk)]=mseqcorrelate_d_s(d-0.001,d+0.001,1,1,Nt,Nr,1,K,Kg,Rx_data,nblk,doppler_scale,seq(:,1).');
%     doppler_scale(nblk)=doppler_scale(nblk)+(0.005-doppler_scale(nblk));
%     doppler_scale(nblk)=d;
    [I,D]=numden(sym(rats(doppler_scale(nblk))));
    I=eval(I);D=eval(D);
    clear Rx_block;
    for i=1:Nr
%         Rx_block(i,:)=interp1((0:length(Rx_blockt(i,:))-1),Rx_blockt(i,:),(0:length(Rx_blockt(i,:))-1)/(1+doppler_scale(nblk)),'spline');
        Rx_block(i,:)=resample(Rx_blockt(i,:),D+I,D);
%         Rx_block(i,:)=Rx_blockt(i,:);
%         for n=1:length(Rx_block(i,:))
%           Rx_block(i,n)=Rx_block(i,n)*(exp(sqrt(-1)*2*pi*doppler_scale(nblk)*(-13000)*(n-1)/9765.625));%48828.125
%         end
    end
%     figure();
%     plot(abs(Rx_block(1,:)));
    
    y=Rx_block(:,length(seq)+Kg+1: length(seq)+Kg+K);
    ola=Rx_block(:,length(seq)+Kg+1+K: length(seq)+Kg+K+L-1);
    y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
    h= FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, 1, K, L, sc_idx, SNR); 
%     figure();
%     plot(abs(h(1,:)));
    
    LLa_cod= log(0.5)*ones(2*Nt,Nbps*K);
    LLR_info= zeros(Nt,Nbit);
    [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_qpsk_fn_28may11(y,LLa_cod,h,K,1,L,Nt,Nr,Nbps,SNRdB);
    S_Est_Iter= S_Est;%((k-1)*Nt+1: k*Nt,:)
    for nt= 1: Nt
        LLe_cod((nt-1)*2+1,:)= randdeintrlv(LLe_cod((nt-1)*2+1,:), 0); %interleaving before delivering to equalizer
        LLe_cod((nt-1)*2+2,:)= randdeintrlv(LLe_cod((nt-1)*2+2,:), 0);
        %perform MAP convolutional decoding with soft information fed back
        [LLR_info(nt,:) LLe_cod((nt-1)*2+1: nt*2,:)] = MAPConvDecoder_R1(LLe_cod((nt-1)*2+1: nt*2,:),Nbit);
        LLe_cod((nt-1)*2+1,:)= randintrlv(LLe_cod((nt-1)*2+1,:), 0); %interleaving before delivering to equalizer
        LLe_cod((nt-1)*2+2,:)= randintrlv(LLe_cod((nt-1)*2+2,:), 0);
        LLa_cod((nt-1)*2+1: nt*2,:)= LLe_cod((nt-1)*2+1: nt*2,:);
    end
    for nt= 1: Nt
        %BER calculation for turbo detection
        Decod_InfoBit= LLR_info(nt,:)<0;
        ErrNum(nt)= sum(block_bit~=Decod_InfoBit);
        ber(nt)= ErrNum(nt)/Nbit;
        
        Dec_CodBit= randdeintrlv(demod_qpsk(S_Est(nt,:)),0);
        Decod_InfoBit1= vitdec(Dec_CodBit,T214,tblen,'cont','hard');
        Decod_InfoBit1= Decod_InfoBit1(tblen+1: Nbit);
        ErrNum1(nt)= sum(block_bit(1: K-tblen)~=Decod_InfoBit1);
        ErrNum2(nt)=sum(code_bitt(nblk,:)~=Dec_CodBit);
        ber1(nt)= ErrNum1(nt)/Nbit;
        ber_rec(nblk)=ber(nt);
        berraw(nt)=ErrNum2(nt)/(2*Nbit);
        ber_recraw(nblk)=berraw(nt);
    end
    RMSEd1=sqrt(mean((doppler_scale-d).^2));
end