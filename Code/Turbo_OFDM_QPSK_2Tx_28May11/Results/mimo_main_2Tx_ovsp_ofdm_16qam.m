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
osf= 1;  %oversampling factor
Nos= Ns/osf;   %(over)samples per symbol
Nbps= 4;          %16QAM
Nt= 2;            %number of transducers

T214 = poly2trellis(4,[17 13]); %For convolutional encoder
tblen=12;               %For convolutional decoder
rate= 1/2;              %convolutional coding rate

K= 1024;                %number of subcarriers
Npilot= 240;            %number of pilot symbols for each OFDM symbol
Nbit= K*Nbps*rate;      %number of information bits in one OFDM symbol
Nrepeat= 8;             %number of OFDM symbols in each packet

SNRdB= 30;              %signal-to-noise ratio in dB
SNR= 10^(SNRdB/10);     %SNR in linear scale
L= 100;                  %length of channel
Ndata= 45000;           %number of samples in K=1024 segment
if(0)
    Chnn_idx=[1:Nr];        %selected channel indices
else
    Chnn_idx= [1:12];        %selected channel indices
end
Nr= length(Chnn_idx);   %# of receive hydrophones

BER= zeros(1,Nt);
BER_All= zeros(Nrepeat,Nt);

%------------load Tx symbols/bits--------------
load Data_mod_ofdm_16qam_1024_2Tx.mat  %load Tx symbols
HQAM_TxSym= Data_mod_ofdm_16qam_1024; 
clear Data_mod_ofdm_16qam_1024; 
load Bitstream_ofdm_16qam_1024_2Tx.mat %load Tx bits
if(0)
load pilot_2Tx_ofdm.mat %pilot symbols
else
load pilot_2Tx_ofdm_up4.mat %upsampled pilot symbols
end

%------------2Tx, K=1024, data structure-------------
%  Tx signal: 2x(lfmb(1k)+gap(0.3k))+lfmb(1k)+gap(0.3k)+m-seq(511)+gap(189)+
%  OPSK blocks(8*1504)+8PSK blocks(8*1504)+16QAM
%  blocks(8*1504)+2x(gap(0.3k)+lfme(1k))

% ---------load Rx signal and resample it-----------
%load RX_2Tx_OFDM_BB_data_3011554F010_C0_S3 %N_shift_symbol= 5;
%load RX_2Tx_OFDM_BB_data_3011554F010_C0_S4 %no detection
%load RX_2Tx_OFDM_BB_data_3011554F010_C1_S3 %bad

%load RX_2Tx_OFDM_BB_data_3011754F010_C0_S3 %good

load RX_2Tx_OFDM_BB_data_3020156F010_C0_S6 %good


counter= 0;
for N_shift_sample= -6: -6
counter= counter+1;

N_shift_sample

for nblk= [1:8]
    
    nblk
    
    Noffset= 1300*Nt+511+189+(16+nblk-1)*1504; %start of the nblk-th 16QAM block (16)
    N_shift_symbol= 6;
    %N_shift_sample= 0;
        
    %--------------initial channel estimation with pilots--------------
    if(1)
        pilot_temp= zeros(Nt,osf*Npilot);
        %channel estimation with pilots (Chu_sequence)
        Gap=Ns*(Noffset+N_shift_symbol)+ N_shift_sample; 
        Y= RX_data(Chnn_idx,Gap: Nos: Gap+osf*K*Nos-1);        
        pilot_temp(:,1: osf: end)= pilot(:,1: Ns: end); 
        if(0)
            h= TimeDomain_MIMO_ChnnEst_fn(Y(:,1:Npilot),pilot_temp,Nt,Nr,L,Npilot,SNR);
        else
            %h= TimeDomain_MIMO_ChnnEst_fn_new(Y(:,1:Npilot),pilot,Nt,Nr,L,Npilot,SNR);
            h= TimeDomain_MIMO_ChnnEst_fn_new(Y(:,1:osf*Npilot+osf*L-1),pilot_temp,Nt,Nr,osf*L,osf*Npilot,SNR);
        end
        %PlotChannel(h,Nt,Nr,osf*L); 
        %size(h)
        %return;        
    elseif(1)
        Kp= 500;
        Gap=Ns*(Noffset+Npilot+120+N_shift_symbol)+ N_shift_sample; 
        Y= RX_data(Chnn_idx,Gap: Nos: Gap+osf*K*Nos-1); 
        if(1)
            pilot= zeros(Nt,osf*Kp);
            temp= zeros(Nt,K);
            for nt= 1: Nt
                temp(nt,:)= ifft(HQAM_TxSym(nt,(nblk-1)*K+1: nblk*K));
            end
            pilot(:,1:osf:osf*Kp)= temp(:,1: Kp);
        else
            temp= zeros(Nt,Ns*K);
            for nt= 1: Nt
                Data_mod= HQAM_TxSym(nt,(nblk-1)*K+1: nblk*K);
                temp(nt,:)= ifft([Data_mod(1:K/2) zeros(1,(Ns-1)*K)...
                    Data_mod(K/2+1:K)],Ns*K)*sqrt(Ns*K);
            end
            pilot(:,1: Nos: osf*Kp*Nos)= temp(:,1: Ns: Kp*Ns);
            %pilot(:,2:2:end)= 0; %why this is important?!
        end
        h= TimeDomain_MIMO_ChnnEst_fn(Y(:,1:osf*Kp),pilot,Nt,Nr,osf*L,osf*Kp,SNR);
        %PlotChannel(h,Nt,Nr,osf*L); 
        %return     
    else
        Kp= 500;
        Gap=Ns*(Noffset+Npilot+120+N_shift_symbol)+ N_shift_sample; 
        Y= RX_data(Chnn_idx,Gap: Nos: Gap+osf*K*Nos-1); 
        for nt= 1: Nt
            temp(nt,:)= ifft(HQAM_TxSym(nt,(nblk-1)*K+1: nblk*K));
        end
        pilot= temp(:, 1: Kp);
        h1= TimeDomain_MIMO_ChnnEst_fn(Y(:,1:osf:osf*Kp),pilot,Nt,Nr,L,Kp,SNR);
        h= zeros(Nr,Nt*L*osf);
        for nr= 1: Nr
            for nt= 1: Nt
                h(nr,(nt-1)*osf*L+1: nt*osf*L) = resample(h1(nr,(nt-1)*L+1: nt*L), osf, 1);
            end
        end
        %PlotChannel(h,Nt,Nr,osf*L);
        %return
    end  
    
    %starting MIMO equalization for OFDM symbol 
    Gap=Ns*(Noffset+Npilot+120+N_shift_symbol)+ N_shift_sample; 
    %Y= RX_data(Chnn_idx,Gap: Ns: Gap+Ndata*Ns-1); Y=[Y zeros(Nr,K1)];
    Yblk= RX_data(Chnn_idx,Gap: Nos: Gap+osf*K*Nos-1); 
    Yblk(:,1:osf*L-1)= Yblk(:,1:osf*L-1)+RX_data(Chnn_idx,Gap+K*Ns: Nos: Gap+K*Ns+(osf*L-1)*Nos-1); %overlap-add
    HQAM_TxSym_blk= HQAM_TxSym(:,(nblk-1)*K+1: nblk*K); %take one block of Tx symbols
    
    x_Est = MIMO_Ovsp_OFDM_MMSE_Det_qpsk(Yblk,HQAM_TxSym_blk,h,K,osf,L,Nt,Nr,Nbps,SNRdB);
    
    if(0)
    %scatterplot(x_Est(1,:)); title('Tx 1');
    %scatterplot(x_Est(2,:)); title('Tx 2');
    scatterplot([x_Est(1,:) x_Est(2,:)]); title('Tx 1 & 2');
    end
    
%     for nt= 1: Nt
%         %perform MAP convolutional decoding with soft information fed back
%         if(0)
%             [LLR_info(nt,:) LLe_cod((nt-1)*2+1: nt*2,:)] = MAPConvDecoder_R1(LLe_cod((nt-1)*2+1: nt*2,:),Nbit);
%         else
%             [LLR_info(nt,:) LLe_cod((nt-1)*2+1: nt*2,:) LL_cod((nt-1)*2+1: nt*2,:)] = MAPConvDecoder_plot(LLe_cod((nt-1)*2+1: nt*2,:),Nbit);
%         end
%     end
% 
    ErrNum1= zeros(Nt,1);
    if(1)
        for nt= 1: Nt
            %BER calculation
            Dec_CodBit= randdeintrlv(demod_mqam(x_Est(nt,:),16),0);
            Decod_InfoBit1= vitdec(Dec_CodBit,T214,tblen,'cont','hard'); 
            Decod_InfoBit1= Decod_InfoBit1(tblen+1: Nbit);
            ErrNum1(nt)= sum(Bitstream_ofdm_16qam_1024(nt,(nblk-1)*K*Nbps/2+1: nblk*K*Nbps/2-tblen)~=Decod_InfoBit1);
            BER(nt)= ErrNum1(nt)/Nbit;
        end
    end
%         
%     %x_Est_all= x_Est;
%     %eval(['save ..\Plotting\x_Est_2Tx_OPSK_K1024_Pkt' int2str(nblk) '_20dB_Iter' int2str(k) '.mat x_Est_all']);
%     eval(['save ..\Plotting\x_Equ_Soft_2Tx_OPSK_K1024_Pkt' int2str(nblk) '_20dB_Iter' int2str(k) '.mat x_Equ_Soft']);        
%     soft_x_MAP= SoftSymCal_OPSK(LL_cod,Nt,K);
%     eval(['save ..\Plotting\soft_x_MAP_2Tx_OPSK_K1024_Pkt' int2str(nblk) '_20dB_Iter' int2str(k) '.mat soft_x_MAP']);
%  
      %BER
      %BER_All(nblk,:)= BER;
      ErrNum1
      ErrNum_All(nblk,counter)= sum(ErrNum1);

end

%BER_All
%ErrNum_All(:,counter)= sum(BER_All*1024,2);
ErrNum_All_Sum= sum(ErrNum_All)

end

BER_mean= ErrNum_All_Sum/(K*Nt*rate*Nbps)/Nrepeat %4 means 16QAM modulation

return




