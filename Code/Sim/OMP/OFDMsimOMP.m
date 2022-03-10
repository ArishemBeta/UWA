clc;
clear all;
close all;

MMode='QPSK';
K=1024*0.5^3;
B=9765.625*0.5^3;
% K=1024;
% B=9765.625;
if(strcmp(MMode,'QPSK'))
    Nbps=2;
elseif(strcmp(MMode,'16QAM'))
    Nbps=4;
end
Nt=1;
Nr=1;
Ns=16;
T214 = poly2trellis(4,[17 13]); %For convolutional encoder
tblen=12;               %For convolutional decoder
rate= 1/2;              %convolutional coding rate
sc_idx= [1:2:K];        %indices of subcarriers for channel estimation
Nbit= K*Nbps*rate;      %number of information bits in one OFDM symbol
SNRdB= 10;              %signal-to-noise ratio in dB
SNR= 10^(SNRdB/10);     %SNR in linear scale
L= 40;                  %length of channel
Kg= 3*L;                %gap between OFDM symbol or between OFDM symbol and pilot block
Nfrm=1;                 % frame
Nblock=10;              % block per frame

%------------generate Tx bits--------------
P_bit=randi([0 1],1,Nbit*Nblock*Nfrm);
% load P_bitQPSK.mat;
bit_seq=reshape(P_bit,Nbit,length(P_bit)/Nbit);

%------------channel coding--------------
for i=1:length(P_bit)/Nbit
code_bitt(i,:)=convenc(bit_seq(:,i),T214);
code_bit(i,:)=randintrlv(code_bitt(i,:),0);
%code_bit=code_bitt;
end


if (strcmp(MMode,'QPSK'))
    data_sym_t=qpsk(code_bit).';
elseif (strcmp(MMode,'16QAM'))
    for i=1:Nblock
        data_sym_t(i,:)=mqam(code_bit(i,:),16);
    end
    data_sym_t=data_sym_t.';
end

%------------频域插零----------
% data_sym=data_sym_t;
data_sym=zeros(K*Ns,Nblock);
for i=1:Nblock
    for j=1:K/2
        data_sym(j,i)=data_sym_t(j,i);
    end
    for j=K/2+1:K
        data_sym(j+(Ns-1)*K,i)=data_sym_t(j,i);
    end
end
% scatterplot(data_sym(:,1));

%------------IFFT--------------
ifft_data=ifft(data_sym);

%------------ZP--------------
Tx_block=[ifft_data;zeros(Kg*Ns,Nfrm*Nblock)];

%------------P/S--------------
Tx_data=reshape(Tx_block,[],1).';
% figure(1);


%------------转通带-----------
% lowpass(Tx_data,B,B*Ns);
for i=1:length(Tx_data)
    Tx_data(i)=real(Tx_data(i)*(exp(sqrt(-1)*2*pi*13000*(i-1)/(B*Ns))));%
end
% plot(Tx_data);

%------------channel--------------
Npath=2;
UWAchannel=UWAchannel_generation(Npath,Nblock*(K+Kg)/B,1/(Ns*B),10000*L/B,1000/B,0.0030*B,B);
Rx_data=zeros(1,Nblock*(K+Kg)*Ns+L*Ns);
for t=1:height(UWAchannel)
    for p=1:Npath
        if(round(t-UWAchannel(t,p*2)*Ns)<=0)
            Rx_data(t)=0;
        else
            Rx_data(t)=Rx_data(t)+UWAchannel(t,p*2-1)*...%round(t-UWAchannel(t,p*2)*Ns)
            Tx_data(round(t-UWAchannel(t,p*2)*Ns));
        end
    end
end
Rx_data=awgn(Rx_data,SNRdB,'measured');
% Rx_data=Tx_data;
% figure();
% plot(Rx_data);

%------------转基带-----------
for i=1:length(Rx_data)
    Rx_data(i)=Rx_data(i)*(exp(sqrt(-1)*2*pi*(-13000)*(i-1)/(B*Ns)));%
end
Rx_data=lowpass(Rx_data,B*0.6,B*Ns);%



%------------S/P--------------
% BER=zeros(1,Nblock);
doppler_scale=zeros(1,Nblock);
cfo=zeros(1,Nblock);
for nblk=1: Nblock
    nblk
    Noffset=(nblk-1)*(K+Kg)-floor(doppler_scale(max(1,nblk-1))*(K+Kg)*(nblk-1)); %
    Gap=Ns*Noffset;
    block_symbol=data_sym_t(:,nblk).';
    pilot_symbol=block_symbol(sc_idx);
    block_bit=bit_seq(:,nblk).';
    pilot_bit=block_bit(sc_idx);
    Rx_block=Rx_data(:,Gap+1:Gap+round((K+Kg)*Ns));%/(1-0.005)                             5~0.0005    30~0.003
%     figure();
%     plot(abs(Rx_block(1,:)));
    
    if (strcmp(MMode,'QPSK'))
        [doppler_scale(nblk),cfo(nblk)]=D2SearchQPSK(-0.0034,-0.0026,0,2,1,1,1,sc_idx,Nt,Nr,Ns,K,L,Rx_block,pilot_symbol,block_symbol,SNR,SNRdB,Nbps,B);
    elseif (strcmp(MMode,'16QAM'))
        [doppler_scale(nblk),cfo(nblk)]=D2SearchQAM(-0.0031,-0.0029,-1,1,1,1,1,sc_idx,Nt,Nr,Ns,K,L,Rx_block,pilot_symbol,block_symbol,SNR,SNRdB,Nbps,B);
    end
%     doppler_scale(nblk)=-0.012;
    for i=1:Nr
        Rx_block(i,:)=interp1((0:length(Rx_block(i,:))-1),Rx_block(i,:),(0:length(Rx_block(i,:))-1)/(1+doppler_scale(nblk)),'spline');
    end
    for i=1:Nr
        for n=1:length(Rx_block(i,:))
            Rx_block(i,n)=Rx_block(i,n)*(exp(sqrt(-1)*2*pi*(cfo(nblk)-13000*doppler_scale(nblk))*(n-1)/(B*Ns)));%48828.125
        end
    end
    
    y=Rx_block(:,1: Ns: Ns*K);
    ola=Rx_block(:,1+K*Ns: Ns: (K+L-1)*Ns);
    y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);

%     z=fft(y).';
%     H=OMP(z,Npath,1,B,Ns,K,sc_idx,block_symbol.',0.00001,0.00001,cfo(nblk));
%     S_Est=(H\z).';

%     if (strcmp(MMode,'QPSK'))
%         Dec_CodBit= randdeintrlv(demod_qpsk(S_Est),0);
%     elseif (strcmp(MMode,'16QAM'))
%         Dec_CodBit= randdeintrlv(demod_mqam(S_Est,16),0);
%     end
%     Decod_InfoBit1= vitdec(Dec_CodBit,T214,tblen,'cont','hard');
%     Decod_InfoBit1= Decod_InfoBit1(tblen+1: Nbit);
%     ErrNum1= sum(block_bit(1: Nbit-tblen)~=Decod_InfoBit1);
%     ErrNum2=sum(code_bitt(nblk,:)~=Dec_CodBit);
%     ber1= ErrNum1/Nbit;
%     ber_rec(nblk)=ber1;
%     berraw=ErrNum2/(2*Nbit);
%     ber_recraw(nblk)=berraw;

    h= FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, Ns, K, L, sc_idx, SNR); 
    figure();
    plot(abs(h(1,:)));
    
    LLa_cod= log(0.5)*ones(2*Nt,Nbps*K);
    LLR_info= zeros(Nt,Nbit);

    if (strcmp(MMode,'QPSK'))
        [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_qpsk_fn_28may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB);    
    elseif (strcmp(MMode,'16QAM'))
        [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_16qam_fn_31may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB);
    end

    scatterplot(S_Est);
    S_Est_Iter= S_Est;%((k-1)*Nt+1: k*Nt,:)

    BER_cost_d=0;
    symbol_est=S_Est;
    pilot_symbol_est=symbol_est(:,sc_idx);
    symbol_err=pilot_symbol-pilot_symbol_est;
    for nt=1:Nt
        for k=1:length(sc_idx)
            BER_cost_d=BER_cost_d+(symbol_err(nt,k)*conj(symbol_err(nt,k)));
        end
    end

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
        if (strcmp(MMode,'QPSK'))
            Dec_CodBit= randdeintrlv(demod_qpsk(S_Est(nt,:)),0);
        elseif (strcmp(MMode,'16QAM'))
            Dec_CodBit= randdeintrlv(demod_mqam(S_Est(nt,:),16),0);
        end
        
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