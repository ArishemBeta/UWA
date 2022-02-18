 %Purpose: Process SPACE08 OFDM data: 2Tx,K=1024,QPSK
%Author: Jun Tao
%
%
clc
clear all
close all
addpath ..\RXSyncData_SPACE08\2Tx_OFDM  %save Sync BB signal
addpath ..\RefSig %LFM, RRC filter etc.
addpath ..\
addpath ..\SPACE08_Moving
Fs=43669.5;  %39.0625 kHz, Sampling Frequency
Fc=13e3;     %carrier frequency  (Hz)
Fb=9765.625;     %9.7656 kHz, Signal Bandwidth
Ns=5;
% Ns=4;
% Ns=1;
Nbps= 2;          %QPSK
Nt= 1;            %number of transducers
T214 = poly2trellis(4,[17 13]); %For convolutional encoder
tblen=12;               %For convolutional decoder
rate= 1/2;              %convolutional coding rate
K= 1024;                %number of subcarriers
Npilot= 90;            %number of pilot symbols for each OFDM symbol
Kg= 120;                %gap between OFDM symbol or between OFDM symbol and pilot block
sc_idx= [1:4:K];        %indices of subcarriers for channel estimation
Nbit= K*Nbps*rate;      %number of information bits in one OFDM symbol
Nrepeat= 12;             %number of OFDM symbols in each packet
SNRdB= 10;              %signal-to-noise ratio in dB
SNR= 10^(SNRdB/10);     %SNR in linear scale
L= 80;                 %length of channel

%Chnn_idx=[1];                              %01 phone
% Chnn_idx=[1 2];                           %02 phones
% Chnn_idx=[1 3 5 7];                      %04 phones
% Chnn_idx=[1 2 3 6 7 8];                  %06 phones
Chnn_idx=[1 2 3 4 5 6 7 8];             %08 phones
Nr= length(Chnn_idx);   %# of receive hydrophones

LL_min= -1e5;
pilot_LLR= 0; %1: using pilot LLR; 0: not using pilot LLR
Enh= 0; %1: enhanced equalization; 0: original equalization
Niter= 1;

%------------load Tx symbols/bits--------------
load One_Tx_Ant1_OFDM_QPSK_K1024.mat  %load Tx symbols
QPSK_TxSym= TxSym_1Tx;
InfoBit=TxInfoBit_1Tx;
CodBit=TxCodIntrlvBit_1Tx;
clear One_Tx_Ant1_OFDM_QPSK_K1024;
load pilot_1Tx_OFDM.mat;
%load Bitstream_ofdm_qpsk_1024_2Tx.mat %load Tx bits

% [packet_name, Sync]=list200;
% packet_error= [];

doppler_scale=zeros(1,12);
cfo=zeros(1,12);
for packet_idx= [1] %[1 2 3 4 5 6 8 11 14 20 23 26 29 32]
    
    BB_name= ['RX_1Tx_OFDM_1024_BB_data_RemusRun01S05F10']
    eval(['load ' BB_name]);
    Dopplera=zeros(1,8);
    for i=1:Nr
        RX_data_n(i,:)=interp1([0:length(RX_data(Chnn_idx(i),:))-1],RX_data(Chnn_idx(i),:),[0:855/956:length(RX_data(Chnn_idx(i),:))-1],'spline');
        RX_data_nn(i,:)=interp1((0:length(RX_data_n(i,:))-1),RX_data_n(i,:),(0:length(RX_data_n(i,:))-1)/(1-0.001312),'spline');
    end
    N_shift_symbol= 10;
    N_shift_sample= 0;
    
    ErrNum= zeros(Nt,Niter);
    ErrNum1= zeros(Nt,Niter);
    BER= zeros(Nt,Niter);
    BER1= zeros(Nt,Niter);
    
    ErrNum_All= zeros(Niter,Nrepeat);
    ErrNum1_All= zeros(Niter,Nrepeat);
    BER_All= zeros(Niter,Nrepeat);
    BER1_All= zeros(Niter,Nrepeat);
    
    ErrNum_All_Tx= zeros(Niter*Nt,Nrepeat);
    
    S_Est_Iter= zeros(Nt*Niter,K);
    
    for nblk= [1: 12]%Nrepeat]
        
        nblk
        
        Noffset= 1300*Nt+511+189+(nblk-1)*1354; %start of the nblk-th block
        Gap=Ns*(Noffset+N_shift_symbol)+ N_shift_sample;
        block_symbol=QPSK_TxSym(nblk,:);
        pilot_symbol=block_symbol(sc_idx);
        block_bit=InfoBit(nblk,:);
        pilot_bit=block_bit(sc_idx);
%         RX_block=RX_data(Chnn_idx,Gap:Gap+K*Ns+120*Ns+Npilot*Ns+Kg*Ns-1);

        %--------------initial channel estimation with pilots--------------
        if(0)
            pilot_temp= zeros(Nt,Npilot);
            %channel estimation with pilots (Chu_sequence)
%             Y= RX_data(Chnn_idx,Gap: Ns: Gap+K*Ns-1);
            RX_block=RX_data_nn(:,Gap:Gap+K*Ns++Npilot*Ns+Kg*Ns+L*Ns-1);
            Y=RX_block(:,1: Ns: K*Ns);
            pilot_temp(:,1: end)= pilot(:,1: end);
            h= TimeDomain_MIMO_ChnnEst_fn_new(Y(:,1:Npilot+L-1),pilot_temp,Nt,Nr,L,Npilot,SNR/L/Nt);
        else
            Gap= Gap+(Npilot+Kg)*Ns; %jump over the pilot block      
            RX_block=RX_data_nn(:,Gap:Gap+K*Ns+L*Ns-1);
%             cfo=[19.1,19.17,19.01,19.05,19.16,19.32,19.26,19.1,19.1,19.16,19,19.1;];
%             for i=1:Nr
%                 for n=1:length(RX_block(i,:))
%                     RX_block(i,n)=RX_block(i,n)*(exp(sqrt(-1)*2*pi*cfo(nblk)*(n-1)/48828.125));
%                 end
%             end
%             doppler_scale(nblk) = DopplerScaleEstimation(1,1,1,1,1,sc_idx,Nt,Nr,Ns,K,L,RX_block,pilot_symbol,block_symbol,pilot_bit,block_bit);
%             doppler_scale(nblk) = DopplerScaleIteration(-0.002,-0.001,1,1,1,sc_idx,Nt,Nr,Ns,K,L,RX_block,pilot_symbol,block_symbol,pilot_bit,block_bit);
%             for i=1:Nr
%                 RX_block(i,:)=interp1((0:length(RX_block(i,:))-1),RX_block(i,:),(0:length(RX_block(i,:))-1)/(1+doppler_scale(nblk)),'spline');
%             end
%             cfo(nblk) = CFOEstimation(1,1,1,1,1,sc_idx,Nt,Nr,Ns,K,L,RX_block,pilot_symbol,block_symbol,pilot_bit,block_bit,pilot);
            cfo(nblk) = CFOIteration(18,22,1,1,1,sc_idx,Nt,Nr,Ns,K,L,RX_block,pilot_symbol,block_symbol,pilot_bit,block_bit);
            for i=1:Nr
                for n=1:length(RX_block(i,:))
                    RX_block(i,n)=RX_block(i,n)*(exp(sqrt(-1)*2*pi*cfo(nblk)*(n-1)/48828.125));
                end
            end
            y=RX_block(:,1: Ns: K*Ns);
            ola=RX_block(:,1+K*Ns: Ns: K*Ns+(L-1)*Ns);
            y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);

%             Gap= Gap+(Npilot+Kg)*Ns; %jump over the pilot block
%             Y= RX_data(Chnn_idx,Gap: Ns: Gap+K*Ns-1); %the k-th received block
%             Y(:,1:L-1)= Y(:,1:L-1)+RX_data(Chnn_idx,Gap+K*Ns: Ns: Gap+K*Ns+(L-1)*Ns-1); %overlap-add
            FDSym= QPSK_TxSym(nblk,:);
            h= FreqDomain_MIMO_ChnnEst_fn_26May11(y, FDSym, Nt, Nr, Ns, K, L, sc_idx, SNR);
        end
        %check CIR estimation
        if(0)
            PlotChannel(h,Nt,Nr,L);
            %save Y.mat Y
            %save H.mat H
            return
        end
        %starting MIMO equalization for OFDM symbol
%         Gap=Ns*(Noffset+Npilot+Kg+N_shift_symbol)+ N_shift_sample;
%         RX_block=RX_block(:,Gap:Gap+K*Ns+L*Ns-1);
%         y=RX_block(:,1051: Ns: K*Ns+1050);
%         ola=RX_block(:,1051+K*Ns: Ns: K*Ns+(L-1)*Ns+1050);
%         y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
        %Y= RX_data(Chnn_idx,Gap: Ns: Gap+Ndata*Ns-1); Y=[Y zeros(Nr,K1)];
%         Yblk= RX_data(Chnn_idx,Gap: Ns: Gap+K*Ns-1);
%         Yblk(:,1:L-1)= Yblk(:,1:L-1)+RX_data(Chnn_idx,Gap+K*Ns: Ns: Gap+K*Ns+(L-1)*Ns-1); %overlap-add
        Yblk=y;
        if(0) %check the spectrum of Rx signal
            figure(1000);
            plot(abs(fft(Yblk(1,:)))); xlim([0 2047]);
            return
        end
        
        %initial LL for coded bits and information bits
        LLa_cod= log(0.5)*ones(2*Nt,Nbps*K);
        LLR_info= zeros(Nt,Nbit);
        
        if(pilot_LLR)
            %replace estimated LLs of pilot bits with their true LLs.
            LLa_p= zeros(2*Nt,Nbps*length(sc_idx));
            for nt= 1: Nt
                pilot_bit= demod_qpsk(QPSK_TxSym(nt,(nblk-1)*K+sc_idx));
                pilot_bit(pilot_bit==0)=0;
                pilot_bit(pilot_bit==1)=LL_min;
                LLa_p((nt-1)*2+1,:)= pilot_bit; %LL for bit equal to 0
                LLa_p((nt-1)*2+2,:)= LL_min-pilot_bit; %LL for bit equal to 1
            end
            %normalization
            for nt= 1: Nt
                for n= 1: length(sc_idx)*Nbps
                    norm1(n)= logsum(LLa_p((nt-1)*2+1: nt*2,n));
                end
                LLa_p((nt-1)*2+1: nt*2,:)= LLa_p((nt-1)*2+1: nt*2,:)-ones(2,1)*norm1;
            end
        end
        
        for k= 1: Niter
            
            k
            
            if(pilot_LLR)
                %replace estimated LLs of pilot bits with their true LLs.
                for nbps= 1: Nbps
                    LLa_cod(:,(sc_idx-1)*Nbps+nbps)= LLa_p(:,nbps:Nbps:end);
                end
            end
            if(~Enh) %original soft equalization
                [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_qpsk_fn_28may11(Yblk,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB);
            else  %enhanced soft equalization
                [S_Est LLe_cod]= MIMO_OFDM_EnhSoftEqu_qpsk_fn_04June11(Yblk,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB);
            end
            S_Est_Iter((k-1)*Nt+1: k*Nt,:)= S_Est;
            if(0)
                scatterplot(S_Est(1,:)); title('Tx 1');
%                 scatterplot(S_Est(2,:)); title('Tx 2');
                %return
            end
            
            if(pilot_LLR)
                %replace estimated LLs of pilot bits with their true LLs.
                for nbps= 1: Nbps
                    LLe_cod(:,(sc_idx-1)*Nbps+nbps)= LLa_p(:,nbps:Nbps:end);
                end
            end
            for nt= 1: Nt
                LLe_cod((nt-1)*2+1,:)= randdeintrlv(LLe_cod((nt-1)*2+1,:), 0); %interleaving before delivering to equalizer
                LLe_cod((nt-1)*2+2,:)= randdeintrlv(LLe_cod((nt-1)*2+2,:), 0);
                
                %perform MAP convolutional decoding with soft information fed back
                [LLR_info(nt,:) LLe_cod((nt-1)*2+1: nt*2,:)] = MAPConvDecoder_R1(LLe_cod((nt-1)*2+1: nt*2,:),Nbit);
                
                %LLRe_cod_temp1(nt,:)= LLRe_cod(nt,:);
                LLe_cod((nt-1)*2+1,:)= randintrlv(LLe_cod((nt-1)*2+1,:), 0); %interleaving before delivering to equalizer
                LLe_cod((nt-1)*2+2,:)= randintrlv(LLe_cod((nt-1)*2+2,:), 0);
                LLa_cod((nt-1)*2+1: nt*2,:)= LLe_cod((nt-1)*2+1: nt*2,:);
            end
            
            if(1)
                for nt= 1: Nt
                    %BER calculation for turbo detection
                    Decod_InfoBit= LLR_info(nt,:)<0;
                    ErrNum(nt,k)= sum(block_bit~=Decod_InfoBit);
                    BER(nt,k)= ErrNum(nt,k)/Nbit;
                    
                    Dec_CodBit= randdeintrlv(demod_qpsk(S_Est(nt,:)),0);
                    Decod_InfoBit1= vitdec(Dec_CodBit,T214,tblen,'cont','hard');
                    Decod_InfoBit1= Decod_InfoBit1(tblen+1: Nbit);
                    ErrNum1(nt,k)= sum(block_bit(1: K-tblen)~=Decod_InfoBit1);
                    BER1(nt,k)= ErrNum1(nt,k)/Nbit;
                end
            end
        end %for loop of k
        ErrNum_All(:,nblk)= sum(ErrNum);
        ErrNum1_All(:,nblk)= sum(ErrNum1);
        for nt= 1: Nt
            ErrNum_All_Tx((nt-1)*Niter+1: nt*Niter,nblk)= ErrNum1(nt,:);
        end
    end
    
    BER_All= ErrNum1_All/Nbit
    ErrNum_All_Tx=[ErrNum_All_Tx sum(ErrNum_All_Tx,2)]
    %name= ['ErrNum_QPSK_K1024_' packet_name(packet_idx,:)];
    %eval(['save Results\' name '.mat ErrNum_All_Tx']);
    %aa= ErrNum_All_Tx(1:Niter,end)+ErrNum_All_Tx(Niter+1: Nt*Niter,end)
    %packet_error= [packet_error aa];
end

return


