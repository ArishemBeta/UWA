%Created by Jun Tao at MST on 27June11
%
clc; clear all; close all;
P_1024= 4; L_1024= 100;  Chnn_idx_1024= 1:11:12;  Dist= 1; SNRdB= 20; Niter= 8;
main_mimo_ofdm_2Tx_qpsk_1024_EnhChnnEst_batch(P_1024,L_1024,Chnn_idx_1024,Dist,SNRdB,Niter);

return

%--------------------------------------------------------------------------
%channel estimaiton: 200m (packet 14), 1000m (packet 11)
clc; clear all; close all;
P_1024= 4; L_1024= 100;  Chnn_idx_1024= 1:12;  Dist= 1; SNRdB= 20; Niter= 8;
main_mimo_ofdm_2Tx_qpsk_1024_EnhChnnEst_batch(P_1024,L_1024,Chnn_idx_1024,Dist,SNRdB,Niter);

clc; clear all; close all;
P_1024= 4; L_1024= 100;  Chnn_idx_1024= 1:12;  Dist= 1; SNRdB= 20; Niter= 8;
main_mimo_ofdm_2Tx_8psk_1024_EnhChnnEst_batch(P_1024,L_1024,Chnn_idx_1024,Dist,SNRdB,Niter); 

clc; clear all; close all;
P_1024= 4; L_1024= 100;  Chnn_idx_1024= 1:12;  Dist= 1; SNRdB= 20; Niter= 8;
main_mimo_ofdm_2Tx_16qam_1024_EnhChnnEst_batch(P_1024,L_1024,Chnn_idx_1024,Dist,SNRdB,Niter); 

%--------------------------------------------------------------------------
clc; clear all; close all;
P_2048= 4; L_2048= 100; Chnn_idx_2048= 1:12; Dist= 1; SNRdB= 20; Niter= 8;
main_mimo_ofdm_2Tx_qpsk_2048_EnhChnnEst_batch(P_2048,L_2048,Chnn_idx_2048,Dist,SNRdB,Niter);

clc; clear all; close all;
P_2048= 4; L_2048= 100; Chnn_idx_2048= 1:12; Dist= 1; SNRdB= 20; Niter= 8;
main_mimo_ofdm_2Tx_8psk_2048_EnhChnnEst_batch(P_2048,L_2048,Chnn_idx_2048,Dist,SNRdB,Niter);

clc; clear all; close all;
P_2048= 4; L_2048= 100; Chnn_idx_2048= 1:12; Dist= 1; SNRdB= 20; Niter= 8;
main_mimo_ofdm_2Tx_16qam_2048_EnhChnnEst_batch(P_2048,L_2048,Chnn_idx_2048,Dist,SNRdB,Niter);

%--------------------------------------------------------------------------
clc; clear all; close all;
P_4096= 8; %pilot spacing
L_4096= 100;  %channel length
Chnn_idx_4096= 1:12;  %channel index
Dist= 1; %1: 200m, 0: 1000m
SNRdB= 20; %SNR in dB
Niter= 8;
main_mimo_ofdm_2Tx_qpsk_4096_EnhChnnEst_batch(P_4096,L_4096,Chnn_idx_4096,Dist,SNRdB,Niter);

clc; clear all; close all;
P_4096= 8; %pilot spacing
L_4096= 100;  %channel length
Chnn_idx_4096= 1:12;  %channel index
Dist= 1; %1: 200m, 0: 1000m
SNRdB= 20; %SNR in dB
Niter= 8;
main_mimo_ofdm_2Tx_8psk_4096_EnhChnnEst_batch(P_4096,L_4096,Chnn_idx_4096,Dist,SNRdB,Niter);

clc; clear all; close all;
P_4096= 8; %pilot spacing
L_4096= 100;  %channel length
Chnn_idx_4096= 1:12;  %channel index
Dist= 1; %1: 200m, 0: 1000m
SNRdB= 20; %SNR in dB
Niter= 8;
main_mimo_ofdm_2Tx_16qam_4096_EnhChnnEst_batch(P_4096,L_4096,Chnn_idx_4096,Dist,SNRdB,Niter);

return

