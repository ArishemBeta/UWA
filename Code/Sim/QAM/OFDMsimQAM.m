clc;
clear all;
close all;

K=1024;
Nbps=4; %QPSK
Nt=1;
Nr=8;
T214 = poly2trellis(4,[17 13]); %For convolutional encoder
tblen=12;               %For convolutional decoder
rate= 1/2;              %convolutional coding rate
sc_idx= [1:4:K];        %indices of subcarriers for channel estimation
Nbit= K*Nbps*rate;      %number of information bits in one OFDM symbol
SNRdB= 10;              %signal-to-noise ratio in dB
SNR= 10^(SNRdB/10);     %SNR in linear scale
Kg= 240;                %gap between OFDM symbol or between OFDM symbol and pilot block
L= 80;                  %length of channel
Nfrm=1;                % frame
Nblock=10;                  % block per frame
M=4;

%------------generate Tx bits--------------
% P_bit=randi([0 1],1,Nbit*Nblock*Nfrm);
load P_bitQAM.mat;
bit_seq=reshape(P_bit,length(P_bit)/Nblock,Nblock);

%------------channel coding--------------
for i=1:Nblock
code_bitt(i,:)=convenc(bit_seq(:,i),T214);
code_bit(i,:)=randintrlv(code_bitt(i,:),0);
end

%------------QAM--------------
for i=1:Nblock
data_sym(i,:)=mqam(code_bit(i,:),16);
end
data_sym=data_sym.';
% scatterplot(data_sym);
%------------S/P--------------
% sym_seq=reshape(data_sym,K,length(data_sym)/K);

%------------IFFT--------------
% ifft_data=ifft(sym_seq); 
ifft_data=ifft(data_sym);

%------------ZP--------------
Tx_block=[ifft_data;zeros(Kg,Nfrm*Nblock)];%zeros(Kg,Nfrm*Nblock);

%------------P/S--------------
Tx_data=reshape(Tx_block,[],1);
% figure(1);
% plot(abs(Tx_data));

%------------Rayleigh channel--------------
SampleRate=9765.625;            %采样频率39062.5
MaximumDopplerShift=0.2;                  %Doppler shift
NumTransmitAntennas=Nt;
NumReceiveAntennas=Nr;
RandomStream='mt19937ar with seed';
PathDelays=[0.001,0.002];          %多径延时向量，s0.0015,0.003
AveragePathGains=[0,-6];             %多径信道增益向量，dB
mimochan = comm.MIMOChannel('SampleRate',SampleRate,'MaximumDopplerShift',MaximumDopplerShift,'PathDelays',PathDelays,...
'AveragePathGains',AveragePathGains,'NumTransmitAntennas',NumTransmitAntennas,'NumReceiveAntennas',NumReceiveAntennas,...
'RandomStream','mt19937ar with seed','Seed',8002,'SpatialCorrelationSpecification','None');%Global stream
Rx_data=mimochan(Tx_data).';
d=0.003;
for i=1:Nr
    Rx_data(i,:)=interp1((0:length(Rx_data(i,:))-1),Rx_data(i,:),(0:length(Rx_data(i,:))-1)*(1+d),'spline');
end
Rx_data=awgn(Rx_data,SNRdB,'measured');
c=-5;% figure(2);
plot(abs(Rx_data(1,:)));
% for i=1:Nr
%     for n=1:length(Rx_data(i,:))
%         Rx_data(i,n)=Rx_data(i,n)*(exp(sqrt(-1)*2*pi*c*(n-1)/9765.625));%48828.125
%     end
% end
% Rx_data=Tx_data;

%------------S/P--------------
% BER=zeros(1,Nblock);
doppler_scale=zeros(1,Nblock);
cfo=zeros(1,Nblock);
for nblk=1: Nblock %
    nblk
%     Noffset=(nblk-1)*(K+Kg)-floor(doppler_scale(max(1,nblk-1))*1264*(nblk-1)); %start of the nblk-th block
%     Gap=Ns*(Noffset+N_shift_symbol);
    block_symbol=data_sym(:,nblk).';
    pilot_symbol=block_symbol(sc_idx);
    block_bit=bit_seq(:,nblk).';
    pilot_bit=block_bit(sc_idx);
%     Rx_block=Rx_data(:,Noffset+1:Noffset+(K+Kg));
%     doppler_scale(nblk)=0.003;
%     cfo(nblk)=0;
%     cfo=CFOIteration(-10,10,1,1,1,sc_idx,Nt,Nr,1,K,L,Rx_block,pilot_symbol,block_symbol);
%     [doppler_scale(nblk),cfo(nblk)]=D2Search(0.002,0.004,-1,1,1,1,1,sc_idx,Nt,Nr,1,K,L,Rx_block,pilot_symbol,block_symbol);
    [doppler_scale(nblk),cfo(nblk),shift(nblk)]=D3Iteration(0.002,0.004,-3,3,1,1,1,sc_idx,Nt,Nr,1,K,L,Rx_data,nblk,doppler_scale,pilot_symbol,block_symbol,SNRdB);
    Noffset=(nblk-1)*(K+Kg)-floor(doppler_scale(max(1,nblk-1))*1264*(nblk-1));
    Rx_block=Rx_data(:,Noffset+1:Noffset+(K+Kg));
    for i=1:Nr
        Rx_block(i,:)=interp1((0:length(Rx_block(i,:))-1),Rx_block(i,:),(0:length(Rx_block(i,:))-1)/(1+doppler_scale(nblk)),'spline');
    end
%       cfo(nblk) = CFOEstimation(-20,20,1,1,1,sc_idx,Nt,Nr,1,K,L,Rx_block,pilot_symbol,block_symbol);
      for i=1:Nr
          for n=1:length(Rx_block(i,:))
              Rx_block(i,n)=Rx_block(i,n)*(exp(sqrt(-1)*2*pi*cfo(nblk)*(n-1)/9765.625));%48828.125
          end
      end
%     figure();
%     plot(abs(Rx_block(1,:)));
    
    y=Rx_block(:,1: K);
    ola=Rx_block(:,1+K: K+L-1);
    y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
    h= FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, 1, K, L, sc_idx, SNR); 
%     figure();
%     plot(abs(h(1,:)));
    
    LLa_cod= log(0.5)*ones(2*Nt,Nbps*K);
    LLR_info= zeros(Nt,Nbit);
    [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_16qam_fn_31may11(y,LLa_cod,h,K,1,L,Nt,Nr,Nbps,SNRdB);
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
        Dec_CodBit= randdeintrlv(demod_mqam(S_Est(nt,:),16),0);
        Decod_InfoBit1= vitdec(Dec_CodBit,T214,tblen,'cont','hard');
        Decod_InfoBit1= Decod_InfoBit1(tblen+1: Nbit);
        ErrNum1(nt)= sum(block_bit(1: Nbit-tblen)~=Decod_InfoBit1);
        ErrNum2(nt)=sum(code_bitt(nblk,:)~=Dec_CodBit);
        ber1(nt)= ErrNum1(nt)/Nbit;
        ber_rec(nblk)=ber(nt);
        berraw(nt)=ErrNum2(nt)/(2*Nbit);
        ber_recraw(nblk)=berraw(nt);
    end
    
%     for i=1:Nr
%         Rx_block(i,:)=interp1((0:length(Rx_block(i,:))-1),Rx_block(i,:),(0:length(Rx_block(i,:))-1)/(1+doppler_scale),'spline');
%     end
%     cfo(nblk) = CFOEstimation(-20,20,1,1,1,sc_idx,Nt,Nr,1,K,L,Rx_block,pilot_symbol,block_symbol);
%     for i=1:Nr
%         for n=1:length(Rx_block(i,:))
%             Rx_block(i,n)=Rx_block(i,n)*(exp(sqrt(-1)*2*pi*cfo(nblk)*(n-1)/9765.625));%48828.125
%         end
%     end
% 
%     y=Rx_block(:,1: K);
%     ola=Rx_block(:,1+K: K+L-1);
%     y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
%     h= FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, 1, K, L, sc_idx, SNR); 
%     figure();
%     plot(abs(h(1,:)));
%     
%     LLa_cod= log(0.5)*ones(2*Nt,Nbps*K);
%     LLR_info= zeros(Nt,Nbit);
%     [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_qpsk_fn_28may11(y,LLa_cod,h,K,1,L,Nt,Nr,Nbps,SNRdB);
%     S_Est_Iter= S_Est;%((k-1)*Nt+1: k*Nt,:)
%     for nt= 1: Nt
%         LLe_cod((nt-1)*2+1,:)= randdeintrlv(LLe_cod((nt-1)*2+1,:), 0); %interleaving before delivering to equalizer
%         LLe_cod((nt-1)*2+2,:)= randdeintrlv(LLe_cod((nt-1)*2+2,:), 0);
%         %perform MAP convolutional decoding with soft information fed back
%         [LLR_info(nt,:) LLe_cod((nt-1)*2+1: nt*2,:)] = MAPConvDecoder_R1(LLe_cod((nt-1)*2+1: nt*2,:),Nbit);
%         LLe_cod((nt-1)*2+1,:)= randintrlv(LLe_cod((nt-1)*2+1,:), 0); %interleaving before delivering to equalizer
%         LLe_cod((nt-1)*2+2,:)= randintrlv(LLe_cod((nt-1)*2+2,:), 0);
%         LLa_cod((nt-1)*2+1: nt*2,:)= LLe_cod((nt-1)*2+1: nt*2,:);
%     end
%     for nt= 1: Nt
%         %BER calculation for turbo detection
%         Decod_InfoBit= LLR_info(nt,:)<0;
%         ErrNum(nt)= sum(block_bit~=Decod_InfoBit);
%         BER(nt)= ErrNum(nt)/Nbit;
%         Dec_CodBit= randdeintrlv(demod_qpsk(S_Est(nt,:)),0);
%         Decod_InfoBit1= vitdec(Dec_CodBit,T214,tblen,'cont','hard');
%         Decod_InfoBit1= Decod_InfoBit1(tblen+1: Nbit);
%         ErrNum1(nt)= sum(block_bit(1: K-tblen)~=Decod_InfoBit1);
%         BER1(nt)= ErrNum1(nt)/Nbit;
%         BER_rec(nblk)=BER(nt);
%     end
end