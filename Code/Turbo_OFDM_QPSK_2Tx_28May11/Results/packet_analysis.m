%Purpose: Process SPACE08 OFDM data: 2Tx,K=1024,QPSK
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
Nbps= 2;          %QPSK
Nt= 2;            %number of transducers

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
L= 80;                 %length of channel
Ndata= 45000;           %number of samples in K=1024 segment
if(0)
    Chnn_idx=[1:Nr];        %selected channel indices
else
    %Chnn_idx=[1 11];                          %02 phones
    %Chnn_idx=[1 3 7 11];                      %04 phones
    %Chnn_idx=[1 3 5 7 9 11];                  %06 phones
    %Chnn_idx=[1 3 5 6 7 9 11 12];             %08 phones
    %Chnn_idx=[1 3 4 5 6 7 9 10 11 12];        %10 phones
    Chnn_idx=[1:12];                          %12 phones
end
Nr= length(Chnn_idx);   %# of receive hydrophones

LL_min= -1e5; 
pilot_LLR= 1; %1: using pilot LLR; 0: not using pilot LLR
Enh= 1; %1: enhanced equalization; 0: original equalization
Niter= 3;

%------------load Tx symbols/bits--------------
load Data_mod_ofdm_qpsk_1024_2Tx.mat  %load Tx symbols
QPSK_TxSym= Data_mod_ofdm_qpsk_1024; 
clear Data_mod_ofdm_qpsk_1024; 
load Bitstream_ofdm_qpsk_1024_2Tx.mat %load Tx bits
if(0)
load pilot_2Tx_ofdm.mat %pilot symbols
else
load pilot_2Tx_ofdm_up4.mat %upsampled pilot symbols
end

[packet_name, Sync]=list200;
packet_error= [];

chnn_idx= 1;
Rx_data_all=[];

if(0)
    for packet_idx= 5%packet_idx= [2 5 8 11 14 17 20 23 26 29 32 35 38 41 44] 

        BB_name= ['RX_2Tx_OFDM_BB_data_' packet_name(packet_idx,:)]
        eval(['load ' BB_name]);

        Rx_data_all= [Rx_data_all RX_data(chnn_idx,:) zeros(1,10^5)];  
    end
end

if(1)
    [packet_name, Sync]=list1000;
    for packet_idx= 10%[1 3 5 7 9 11 12 13 14 15] %[1 2 4 7 10 11 12 13 14 15]

        BB_name= ['RX_2Tx_OFDM_BB_data_' packet_name(packet_idx,:)]
        eval(['load ' BB_name]);

        Rx_data_all= [Rx_data_all RX_data(chnn_idx,:) zeros(1,10^5)];  
    end
end

figure(1);
plot(abs(Rx_data_all));

%specgram(Rx_data_all,256,Fs);

return


