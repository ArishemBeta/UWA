%Purpose: Process SPACE08 OFDM data: 2Tx,K=1024,OPSK
%Author: Jun Tao
%Place: MST
%Date: June 28, 2010
%Date: July 01, 2010 (oversampled ofdm)
%
%
clc
clear all
close all
  
addpath ..\RXData_Nov1908_Jian %contain original one-minute Rx PB signal
addpath ..\RXSyncData_SPACE08\2Tx_OFDM  %save Sync BB signal
addpath ..\RXPackage_SPACE08\2Tx_OFDM   %contain divided PB signal (four packages)
addpath ..\RefSig %LFM, RRC filter etc.

Fs=1e7/256;  %39.0625 kHz, Sampling Frequency
Fc=13e3;     %carrier frequency  (Hz)
Fb=Fs/4;     %9.7656 kHz, Signal Bandwidth
Ns=round(Fs/Fb);  %# of samples per symbol
Nbps= 3;          %OPSK
Nt= 2;            %number of transducers
Nr= 12;           %number of hydrophones

T214 = poly2trellis(4,[17 13]); %For convolutional encoder
tblen=12;               %For convolutional decoder
rate= 1/2;              %convolutional coding rate

K= 1024;                %number of subcarriers
Npilot= 240;            %number of pilot symbols for each OFDM symbol
Kg= 120;                %gap between OFDM symbol or between OFDM symbol and pilot block
sc_idx= [1:3:K];        %indices of subcarriers for channel estimation
Nbit= K*Nbps*rate;      %number of information bits in one OFDM symbol
Nrepeat= 8;             %number of OFDM symbols in each packet

SNRdB= 20;              %signal-to-noise ratio in dB
SNR= 10^(SNRdB/10);     %SNR in linear scale
L= 80;                  %length of channel
Ndata= 45000;           %number of samples in K=1024 segment
if(1)
    Chnn_idx=[1:Nr];        %selected channel indices
else
    %Chnn_idx= [1 3 4 5 6 7 9 10 11 12];    %selected channel indices
    Chnn_idx=[1:12];                          %12 phones
end
Nr= length(Chnn_idx);   %# of receive hydrophones

BER= zeros(1,Nt);
BER_All= zeros(Nrepeat,Nt);

LL_min= -100000; 
pilot_LLR= 1; %1: using pilot LLR; 0: not using pilot LLR
Enh= 1; %1: enhanced equalization; 0: original equalization
Niter= 3;

%------------load Tx symbols/bits--------------
load Data_mod_ofdm_8psk_1024_2Tx.mat  %load Tx symbols
OPSK_TxSym= Data_mod_ofdm_8psk_1024; 
clear Data_mod_ofdm_8psk_1024; 
load Bitstream_ofdm_8psk_1024_2Tx.mat %load Tx bits
if(0)
load pilot_2Tx_ofdm.mat %pilot symbols
else
load pilot_2Tx_ofdm_up4.mat %upsampled pilot symbols
end

if(1)
[packet_name, Sync]=list200;   %[1 2 3 4 5 6 8 11 14 20 23 26 29 32]
else
[packet_name, Sync]=list1000;  %[11]
end
packet_error= [];

for packet_idx= [1 2 3 4 5 6 8 11 14 20 23 26 29 32]
    
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

    for nblk= [1:8]
    
        nblk

        Noffset= 1300*Nt+511+189+(8+nblk-1)*1504; %start of the nblk-th 8PSK block (8)
        Gap=Ns*(Noffset+N_shift_symbol)+ N_shift_sample; 

        %--------------initial channel estimation with pilots--------------
        if(0)
            pilot_temp= zeros(Nt,Npilot);
            %channel estimation with pilots (Chu_sequence)
            Y= RX_data(Chnn_idx,Gap: Ns: Gap+K*Ns-1);        
            pilot_temp(:,1: end)= pilot(:,1: Ns: end); 
            h= TimeDomain_MIMO_ChnnEst_fn_new(Y(:,1:Npilot+L-1),pilot_temp,Nt,Nr,L,Npilot,SNR/L/Nt);   
        else
            Gap= Gap+(Npilot+Kg)*Ns; %jump over the pilot block
            Y= RX_data(Chnn_idx,Gap: Ns: Gap+K*Ns-1); %the k-th received block
            Y(:,1:L-1)= Y(:,1:L-1)+RX_data(Chnn_idx,Gap+K*Ns: Ns: Gap+K*Ns+(L-1)*Ns-1); %overlap-add
            FDSym= OPSK_TxSym(:,(nblk-1)*K+1: nblk*K);
            h= FreqDomain_MIMO_ChnnEst_fn_26May11(Y, FDSym, Nt, Nr, Ns, K, L, sc_idx, SNR);
        end 

        %check CIR estimation
        if(0) 
            PlotChannel(h,Nt,Nr,L);          
            %save Y.mat Y
            %save H.mat H
            return
        end
        
        %starting MIMO equalization for OFDM symbol 
        Gap=Ns*(Noffset+Npilot+Kg+N_shift_symbol)+ N_shift_sample; 
        %Y= RX_data(Chnn_idx,Gap: Ns: Gap+Ndata*Ns-1); Y=[Y zeros(Nr,K1)];
        Yblk= RX_data(Chnn_idx,Gap: Ns: Gap+K*Ns-1); 
        Yblk(:,1:L-1)= Yblk(:,1:L-1)+RX_data(Chnn_idx,Gap+K*Ns: Ns: Gap+K*Ns+(L-1)*Ns-1); %overlap-add
        
        %initial LL for coded bits and information bits
        LLa_cod= log(0.5)*ones(2*Nt,Nbps*K);  
        LLR_info= zeros(Nt,Nbit); 

        if(pilot_LLR)
            %replace estimated LLs of pilot bits with their true LLs.       
            LLa_p= zeros(2*Nt,Nbps*length(sc_idx));
            for nt= 1: Nt
                pilot_bit= demod_8psk(OPSK_TxSym(nt,(nblk-1)*K+sc_idx));
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
            if(~Enh)
                [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_8psk_fn_28may11(Yblk,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB);
            else
                [S_Est LLe_cod]= MIMO_OFDM_EnhSoftEqu_8psk_fn_04June11(Yblk,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB);
            end
            S_Est_Iter((k-1)*Nt+1: k*Nt,:)= S_Est;
            if(0)
            scatterplot(S_Est(1,:))
            scatterplot(S_Est(2,:))
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

            if(pilot_LLR)
                %replace estimated LLs of pilot bits with their true LLs.
                for nbps= 1: Nbps
                    LLe_cod(:,(sc_idx-1)*Nbps+nbps)= LLa_p(:,nbps:Nbps:end);
                end
            end
            for nt= 1: Nt
                %BER calculation for turbo detection
                Decod_InfoBit= LLR_info(nt,:)<0;
                ErrNum(nt,k)= sum(Bitstream_ofdm_8psk_1024(nt,(nblk-1)*Nbit+1: nblk*Nbit)~=Decod_InfoBit);
                BER(nt,k)= ErrNum(nt,k)/Nbit;

                Dec_CodBit= randdeintrlv(demod_8psk(S_Est(nt,:)),0);
                Decod_InfoBit1= vitdec(Dec_CodBit,T214,tblen,'cont','hard'); 
                Decod_InfoBit1= Decod_InfoBit1(tblen+1: Nbit);
                ErrNum1(nt,k)= sum(Bitstream_ofdm_8psk_1024(nt,(nblk-1)*Nbit+1: nblk*Nbit-tblen)~=Decod_InfoBit1);
                BER1(nt,k)= ErrNum1(nt,k)/Nbit;
            end
        end %for loop of k    
        ErrNum_All(:,nblk)= sum(ErrNum);
        ErrNum1_All(:,nblk)= sum(ErrNum1);
        for nt= 1: Nt
            ErrNum_All_Tx((nt-1)*Niter+1: nt*Niter,nblk)= ErrNum(nt,:); 
        end
    end
    BER_All= ErrNum_All/Nbit;
    ErrNum_All_Tx=[ErrNum_All_Tx sum(ErrNum_All_Tx,2)]
    name= ['ErrNum_8PSK_K1024_' packet_name(packet_idx,:)];
    eval(['save Results\' name '.mat ErrNum_All_Tx']);
    aa= ErrNum_All_Tx(1:Niter,end)+ErrNum_All_Tx(Niter+1: Nt*Niter,end)
    packet_error= [packet_error aa];
end

return




