 %Purpose: Process SPACE08 OFDM data: 2Tx,K=1024,QPSK
%Author: Jun Tao
%
%
clc
clear all
close all

addpath DopplerScaleEstimation/Code/RXSyncData_SPACE08/2Tx_OFDM  %save Sync BB signal
addpath DopplerScaleEstimation/Code/RefSig %LFM, RRC filter etc.
addpath DopplerScaleEstimation/Code/Turbo_OFDM_QPSK_2Tx_28May11

Fs=1e7/256;  %39.0625 kHz, Sampling Frequency
Fc=13e3;     %carrier frequency  (Hz)
Fb=Fs/4;     %9.7656 kHz, Signal Bandwidth
Ns=round(Fs/Fb);  %# of samples per symbol
Nbps= 2;          %QPSK
Nt= 2;            %number of transducers

T214 = poly2trellis(4,[17 13]); %For convolutional encoder
tblen=12;               %For convolutional decoder
rate= 1/2;              %convolutional coding rate

K= 1024;                %number of subcarriers
Npilot= 240;            %number of pilot symbols for each OFDM symbol
Kg= 120;                %gap between OFDM symbol or between OFDM symbol and pilot block
sc_idx= [1:4:K];        %indices of subcarriers for channel estimation
Nbit= K*Nbps*rate;      %number of information bits in one OFDM symbol
Nrepeat= 8;             %number of OFDM symbols in each packet

SNRdB= 15;              %signal-to-noise ratio in dB
SNR= 10^(SNRdB/10);     %SNR in linear scale
L= 100;                 %length of channel
%Ndata= 45000;           %number of samples in K=1024 segment

% Chnn_idx=[1];                              %01 phone
% Chnn_idx=[1 11];                           %02 phones
% Chnn_idx=[1 3 7 11];                      %04 phones
% Chnn_idx=[1 3 5 7 9 11];                  %06 phones
% Chnn_idx=[1 3 5 6 7 9 11 12];             %08 phones
% Chnn_idx=[1 3 4 5 6 7 9 10 11 12];        %10 phones
Chnn_idx=[1:12];                          %12 phones
Nr= length(Chnn_idx);   %# of receive hydrophones

LL_min= -1e5;
pilot_LLR= 0; %1: using pilot LLR; 0: not using pilot LLR
Enh= 0; %1: enhanced equalization; 0: original equalization
Niter= 1;

%------------load Tx symbols/bits--------------
load Data_mod_ofdm_qpsk_1024_2Tx.mat  %load Tx symbols
QPSK_TxSym= Data_mod_ofdm_qpsk_1024;
clear Data_mod_ofdm_qpsk_1024;
load Bitstream_ofdm_qpsk_1024_2Tx.mat %load Tx bits

[packet_name, Sync]=list200;
packet_error= [];

for packet_idx= [1] %[1 2 3 4 5 6 8 11 14 20 23 26 29 32]
    
    BB_name= ['RX_2Tx_OFDM_BB_data_' packet_name(packet_idx,:)]
    eval(['load ' BB_name]);
    
    N_shift_symbol= Sync(packet_idx,1);
    N_shift_sample= Sync(packet_idx,2);
    
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
    
    for nblk= [1: 1]%Nrepeat]
        
        nblk
        
        Noffset= 1300*Nt+511+189+(nblk-1)*1504; %start of the nblk-th block
        Gap=Ns*(Noffset+N_shift_symbol)+ N_shift_sample;
        block_symbol=QPSK_TxSym(:,(nblk-1)*K+1: nblk*K);
        pilot_symbol=block_symbol(:,sc_idx);
        block_bit=Bitstream_ofdm_qpsk_1024(:,(nblk-1)*K+1: nblk*K);
        pilot_bit=block_bit(:,sc_idx);

        %         RX_block=RX_data(Chnn_idx,Gap+1:Gap+K*Ns+Kg*Ns+Npilot*Ns+Kg*Ns);

        Gap= Gap+(Npilot+Kg)*Ns; %jump over the pilot block
        RX_block=RX_data(Chnn_idx,Gap+1:Gap+K*Ns+Kg*Ns);

        %------------OMP-------------
        ifOMP=1;
        if(ifOMP)
            yO=RX_block(:,1: Ns: Ns*K).';
            olaO=RX_block(:,1+K*Ns: Ns: (K+L-1)*Ns).';
            yO(1:L-1,:)=(yO(1:L-1,:)+olaO(1:L-1,:));

            z=(fft(yO)./sqrt(height(yO)));
            D=4;
            tic
            H=OMP_test(z,0,Fb,K,sc_idx,block_symbol.',-0.00005,0.00005,0.00001,-1,L,D,SNR);
            toc
            %-------------均衡--------------
            %  S_EstO=((H'*H+N0*eye(K))\H'*z).';
            P=sum(sum(abs(RX_block).^2))/((K+Kg)*Nr*Ns);
            N0=P/(1+SNR);
            tic
%             S_EstO=MIMO_LMMSE_Equalization(reshape(z,1,K*Nr),H,K,[1,1,1,1],Nr,Nt,'QPSK',0,N0);
            S_EstO=MIMO_LMMSE_Equalization(z,H,K,[1,1,1,1],Nr,Nt,'QPSK',0,N0,D);
            toc
            SO=S_EstO;
            for  nt=1:Nt
                Dec_CodBitO= randdeintrlv(demod_qpsk(SO(nt,:)),0);
                Decod_InfoBitO= vitdec(Dec_CodBitO,T214,tblen,'cont','hard');
                Decod_InfoBitO= Decod_InfoBitO(tblen+1: Nbit);
                ErrNum1= sum(Bitstream_ofdm_qpsk_1024(nt,(nblk-1)*K+1: nblk*K-tblen)~=Decod_InfoBitO);
%             ErrNum2=sum(CodBit(nblk,:)~=demod_qpsk(SO));
                ber_recO(nt,nblk)= ErrNum1/(Nbit);
            end
%             ber_recrawO(nblk)=ErrNum2/(Nbit/rate);
            scatterplot(SO(1,:));
            title('OMP',FontSize=20);

            BER_costO=0;
            symbol_errO=block_symbol-SO;
            for nt=1:Nt
                for k=1:width(symbol_errO)
                    BER_costO=BER_costO+(symbol_errO(nt,k)*conj(symbol_errO(nt,k)));
                end
            end
        end


        
%-------------LS-------------
            

%             y=RX_block(:,1441: Ns: K*Ns+1440);
%             ola=RX_block(:,1441+K*Ns: Ns: K*Ns+(L-1)*Ns+1440);
            y=RX_block(:,1: Ns: K*Ns);
            ola=RX_block(:,K*Ns+1: Ns: K*Ns+(L-1)*Ns);
            y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
            
%             Y= RX_data(Chnn_idx,Gap: Ns: Gap+K*Ns-1); %the k-th received block
%             Y(:,1:L-1)= Y(:,1:L-1)+RX_data(Chnn_idx,Gap+K*Ns: Ns: Gap+K*Ns+(L-1)*Ns-1); %overlap-add
            FDSym= QPSK_TxSym(:,(nblk-1)*K+1: nblk*K);
            h= FreqDomain_MIMO_ChnnEst_fn_26May11(y, FDSym, Nt, Nr, Ns, K, L, sc_idx, SNR);
        
        %check CIR estimation
        if(0)
            PlotChannel(h,Nt,Nr,L);
            %save Y.mat Y
            %save H.mat H
            return
        end
        
        %starting MIMO equalization for OFDM symbol
        Gap=Ns*(Noffset+Npilot+Kg+N_shift_symbol)+ N_shift_sample;
        Yblk= RX_data(Chnn_idx,Gap: Ns: Gap+K*Ns-1);
        Yblk(:,1:L-1)= Yblk(:,1:L-1)+RX_data(Chnn_idx,Gap+K*Ns: Ns: Gap+K*Ns+(L-1)*Ns-1); %overlap-add
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
            scatterplot(S_Est(1,:));
            title('LS',FontSize=20);
            
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
                    ErrNum(nt,k)= sum(Bitstream_ofdm_qpsk_1024(nt,(nblk-1)*K+1: nblk*K)~=Decod_InfoBit);
                    BER(nt,k)= ErrNum(nt,k)/Nbit;
                    
                    Dec_CodBit= randdeintrlv(demod_qpsk(S_Est(nt,:)),0);
                    Decod_InfoBit1= vitdec(Dec_CodBit,T214,tblen,'cont','hard');
                    Decod_InfoBit1= Decod_InfoBit1(tblen+1: Nbit);
                    ErrNum1(nt,k)= sum(Bitstream_ofdm_qpsk_1024(nt,(nblk-1)*K+1: nblk*K-tblen)~=Decod_InfoBit1);
                    BER1(nt,k)= ErrNum1(nt,k)/Nbit;
                end
                BER_cost=0;
                symbol_err=block_symbol-S_Est;
                for nt=1:Nt
                    for k=1:length(symbol_err)
                        BER_cost=BER_cost+(symbol_err(nt,k)*conj(symbol_err(nt,k)));
                    end
                end
            end
        end %for loop of k
        
        ErrNum_All(:,nblk)= sum(ErrNum);
        ErrNum1_All(:,nblk)= sum(ErrNum1);
        for nt= 1: Nt
            ErrNum_All_Tx((nt-1)*Niter+1: nt*Niter,nblk)= ErrNum(nt,:);
        end
    end
    
    BER_All= ErrNum_All/Nbit
    ErrNum_All_Tx=[ErrNum_All_Tx sum(ErrNum_All_Tx,2)]
    %name= ['ErrNum_QPSK_K1024_' packet_name(packet_idx,:)];
    %eval(['save Results\' name '.mat ErrNum_All_Tx']);
    %aa= ErrNum_All_Tx(1:Niter,end)+ErrNum_All_Tx(Niter+1: Nt*Niter,end)
    %packet_error= [packet_error aa];
end

return


