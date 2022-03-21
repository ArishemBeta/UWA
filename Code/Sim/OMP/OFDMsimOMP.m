clc;
clear all;
close all;

% parpool('local',16);

MMode='QPSK';
K=1024*0.5^2;
B=9765.625*0.5^0;
ifOMP=1;
% K=1024;
% B=9765.625;

if(strcmp(MMode,'QPSK'))
    Nbps=2;
elseif(strcmp(MMode,'16QAM'))
    Nbps=4;
end
sc=[0,1,0,0,1,1,1,0];
nKr=length(find(sc==0))/length(sc);

Nt=1;
Nr=1;
Ns=16;
T214 = poly2trellis(4,[17 13]);
tblen=12;               
rate= 1/2;              
sc_idx= [2:8:K];                            % pilot index
Nbit_s= K*Nbps*rate;                        % bits in an OFDM symbol including null subcarriers
Nbit_d=Nbit_s*(1-nKr);                      % useful bits in an OFDM symbol
SNRdB= 20;              
SNR= 10^(SNRdB/10);
Lms=8;                                      % length of channel (ms)
L= min(fix(Lms*B/1000),length(sc_idx));     % length of channel
% L=16;
Kg= 3*L;                                    % guard gap
Nfrm=1;                                     % frame
Nblock=10;                                  % block per frame

%------------generate Tx bits--------------
P_bit=randi([0 1],1,Nbit_d*Nblock*Nfrm);
% load P_bitQPSK.mat;
bit_seq=reshape(P_bit,Nbit_d,length(P_bit)/(Nbit_d));

%------------channel coding--------------
for i=1:length(P_bit)/(Nbit_d)
code_bitt(i,:)=convenc(bit_seq(:,i),T214);
code_bit(i,:)=randintrlv(code_bitt(i,:),0);
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
% data_sym=zeros(K*Ns,Nblock);
% for i=1:Nblock
%     for j=1:K/2
%         data_sym(j,i)=data_sym_t(j,i);
%     end
%     for j=K/2+1:K
%         data_sym(j+(Ns-1)*K,i)=data_sym_t(j,i);
%     end
% end
[Tx_sym_ins, Tx_sym]=Txsymbol_arrangement(data_sym_t,K,sc,Ns);
% scatterplot(data_sym(:,1));

%------------IFFT--------------
ifft_data=ifft(Tx_sym_ins);

%------------ZP--------------
Tx_block=[ifft_data;zeros(Kg*Ns,Nfrm*Nblock)];

%------------P/S--------------
Tx_data=reshape(Tx_block,[],1).';

%------------转通带-----------
for i=1:length(Tx_data)
    Tx_data(i)=real(Tx_data(i)*(exp(sqrt(-1)*2*pi*13000*(i-1)/(B*Ns))));%
end

%------------channel--------------
Npath=2;
delay=[1,3,5];
% ddd=-0.00300;
doppler=[-0.00200,-0.00300,-0.002500];
len=8;                                                  %ms
res=0.001;                                               %ms
UWAchannel=UWAchannel_generation(Npath,Nblock*(K+Kg)/B,1/(Ns*B),len,res,delay,doppler);
Rx_data=zeros(1,Nblock*(K+Kg)*Ns+L*Ns);
for t=1:height(UWAchannel)
    for p=1:Npath
        if(round(t-UWAchannel(t,p*2)*Ns*B)<=0)
            Rx_data(t)=0;
        else
            Rx_data(t)=Rx_data(t)+UWAchannel(t,p*2-1)*...%round(t-UWAchannel(t,p*2)*Ns)
            Tx_data(round(t-UWAchannel(t,p*2)*Ns*B));
        end
    end
end
Rx_data=awgn(Rx_data,SNRdB,'measured');

%------------转基带-----------
for i=1:length(Rx_data)
    Rx_data(i)=Rx_data(i)*(exp(sqrt(-1)*2*pi*(-13000)*(i-1)/(B*Ns)));%
end
Rx_data=lowpass(Rx_data,B*0.5,B*Ns);

%------------S/P--------------
doppler_scale=zeros(1,Nblock);
cfo=zeros(1,Nblock);
for nblk=1: 1
    nblk
    Noffset=(nblk-1)*(K+Kg)-floor(doppler_scale(max(1,nblk-1))*(K+Kg)*(nblk-1)); %
    Gap=Ns*Noffset;
    block_symbol=Tx_sym(:,nblk).';
    pilot_symbol=block_symbol(sc_idx);
    block_bit=bit_seq(:,nblk).';
    Rx_block=Rx_data(:,Gap+1:Gap+round((K+Kg)*Ns));%/(1-0.005)                             5~0.0005    30~0.003
    
    doppler_scale(nblk)=mean(doppler(1:Npath));
    cfo(nblk)=0.01;
    if (strcmp(MMode,'QPSK'))
%         [doppler_scale(nblk),cfo(nblk)]=D2SearchQPSK(doppler(1)-0.001,doppler(1)+0.001,-0,0,1,1,1,sc_idx,Nt,Nr,Ns,K,L,Rx_block,pilot_symbol,block_symbol,SNR,SNRdB,Nbps,B);
    elseif (strcmp(MMode,'16QAM'))
        [doppler_scale(nblk),cfo(nblk)]=D2SearchQAM(-0.0031,-0.0029,0,0,1,1,1,sc_idx,Nt,Nr,Ns,K,L,Rx_block,pilot_symbol,block_symbol,SNR,SNRdB,Nbps,B);
    end
    for i=1:Nr
        Rx_block(i,:)=interp1((0:length(Rx_block(i,:))-1),Rx_block(i,:),(0:length(Rx_block(i,:))-1)/(1+doppler_scale(nblk)),'spline');
    end
    for i=1:Nr
        for n=1:length(Rx_block(i,:))
            Rx_block(i,n)=Rx_block(i,n)*(exp(sqrt(-1)*2*pi*(-cfo(nblk)-13000*doppler_scale(nblk))*(n-1)/(B*Ns)));%48828.125
        end
    end

%------------OMP------------
if(ifOMP)
    yO=Rx_block(:,1: Ns: Ns*K);
%     ola=Rx_block(:,1+K*Ns: Ns: (K+L-1)*Ns);
%     y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);

    z=fft(yO).';
    tic
    H=OMP(z,Npath,2,B,Ns,K,sc_idx,block_symbol.',-0.0005,0.0005,0.0005,cfo(nblk),L);
    toc
    S_EstO=(H\z).';
    SO=(diag(repmat(sc,1,K/length(sc)))*S_EstO.').';
    SO(SO==0)=[];

    if (strcmp(MMode,'QPSK'))
        Dec_CodBitO= randdeintrlv(demod_qpsk(SO),0);
    elseif (strcmp(MMode,'16QAM'))
        Dec_CodBitO= randdeintrlv(demod_mqam(SO,16),0);
    end
    Decod_InfoBitO= vitdec(Dec_CodBitO,T214,tblen,'cont','hard');
    Decod_InfoBitO= Decod_InfoBitO(tblen+1: Nbit_d);
    ErrNum1= sum(block_bit(1: Nbit_d-tblen)~=Decod_InfoBitO);
    ErrNum2=sum(code_bitt(nblk,:)~=Dec_CodBitO);
    ber_recO(nblk)= ErrNum1/(Nbit_d);
    ber_recrawO(nblk)=ErrNum2/(Nbit_d/rate);
    scatterplot(SO);

    BER_costO=0;
    pilot_symbol_estO=S_EstO(:,sc_idx);
    symbol_errO=pilot_symbol-pilot_symbol_estO;
    for nt=1:Nt
        for k=1:length(sc_idx)
            BER_costO=BER_costO+(symbol_errO(nt,k)*conj(symbol_errO(nt,k)));
        end
    end
end

%-----------原来的方法-------------
    y=Rx_block(:,1: Ns: Ns*K);
    ola=Rx_block(:,1+K*Ns: Ns: (K+L-1)*Ns);
    y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
    h= FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, Ns, K, L, sc_idx, SNR);
    figure();
    plot(abs(h(1,:)));

    LLa_cod= log(0.5)*ones(2*Nt,Nbps*K);
    LLR_info= zeros(Nt,Nbit_s);

    if (strcmp(MMode,'QPSK'))
        [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_qpsk_fn_28may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB);
    elseif (strcmp(MMode,'16QAM'))
        [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_16qam_fn_31may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB);
    end
    S_Est_Iter= S_Est;%((k-1)*Nt+1: k*Nt,:)

    BER_cost=0;
    pilot_symbol_est=S_Est(:,sc_idx);
    symbol_err=pilot_symbol-pilot_symbol_est;
    for nt=1:Nt
        for k=1:length(sc_idx)
            BER_cost=BER_cost+(symbol_err(nt,k)*conj(symbol_err(nt,k)));
        end
    end

%     for nt= 1: Nt
%         LLe_cod((nt-1)*2+1,:)= randdeintrlv(LLe_cod((nt-1)*2+1,:), 0); %interleaving before delivering to equalizer
%         LLe_cod((nt-1)*2+2,:)= randdeintrlv(LLe_cod((nt-1)*2+2,:), 0);
%         %perform MAP convolutional decoding with soft information fed back
%         [LLR_info(nt,:) LLe_cod((nt-1)*2+1: nt*2,:)] = MAPConvDecoder_R1(LLe_cod((nt-1)*2+1: nt*2,:),Nbit_s);
%         LLe_cod((nt-1)*2+1,:)= randintrlv(LLe_cod((nt-1)*2+1,:), 0); %interleaving before delivering to equalizer
%         LLe_cod((nt-1)*2+2,:)= randintrlv(LLe_cod((nt-1)*2+2,:), 0);
%         LLa_cod((nt-1)*2+1: nt*2,:)= LLe_cod((nt-1)*2+1: nt*2,:);
%     end
    for nt= 1: Nt
%         Decod_InfoBit= LLR_info(nt,:)<0;
%         ErrNum(nt)= sum(block_bit~=Decod_InfoBit);
%         ber(nt)= ErrNum(nt)/Nbit_s;

        S=(diag(repmat(sc,1,K/length(sc)))*S_Est(nt,:).').';
        S(find(S==0))=[];
        scatterplot(S);
        if (strcmp(MMode,'QPSK'))
            Dec_CodBit= randdeintrlv(demod_qpsk(S),0);
        elseif (strcmp(MMode,'16QAM'))
            Dec_CodBit= randdeintrlv(demod_mqam(S,16),0);
        end

        Decod_InfoBit1= vitdec(Dec_CodBit,T214,tblen,'cont','hard');
        Decod_InfoBit1= Decod_InfoBit1(tblen+1: Nbit_d);
        ErrNum1(nt)= sum(block_bit(1: Nbit_d-tblen)~=Decod_InfoBit1);
        ErrNum2(nt)=sum(code_bitt(nblk,:)~=Dec_CodBit);
        ber1(nt)= ErrNum1(nt)/(Nbit_d);
        ber_rec(nblk)=ber1(nt);
        ber_recraw(nblk)=ErrNum2(nt)/(Nbit_d/rate);
    end
    
end
% delete(gcp('nocreate'));