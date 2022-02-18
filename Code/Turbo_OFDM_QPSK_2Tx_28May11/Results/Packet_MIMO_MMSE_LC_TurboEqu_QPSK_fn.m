%-------------------------------------------------
%do block-based turbo equalization for a whole packet
%-------Created by J. Tao at MST on Feb.20,2009----
%
function [x_Est Psi LLRe_cod]=Packet_MIMO_MMSE_LC_TurboEqu_QPSK_fn(Y,QPSK_TxSym,LLRa_cod,Nt,Nr,L,Kp1,Kp2,Nblk,Ngsize,SNRdB);
%
Ndata= 30000;     %number of transmitted data symbols
pcyc= 30;         %less than 20% training
Nbps= 2;          %QPSK

%L= 10;           %length of channel
K2=L;   %Causal part of Feedforward filter (K2 >= 0)
K1=L;   %Anti-causal part of Feedforward filter
K3=0;    %Feedback part, K3 <= K2+L-1
KK= K1+K2+L;

x_Est= QPSK_TxSym(:,1:Kp1); %store pilot symbols
Psi= [];
LLRe_cod= zeros(Nt,Kp1*Nbps); %LLR for pilot symbols(bits)
    
%------------ for first block equalization ---------------------------
%----Initial channel estimation using training symbols at the front of data packet 
SNR= 10^(SNRdB/10);  %SNR in linear scale
h= TimeDomain_MIMO_ChnnEst_fn(Y(:,1:Kp1),QPSK_TxSym(:,1:Kp1),Nt,Nr,L,Kp1,SNR);
%size(h)
%check CIR estimation
if(0) 
    PlotChannel(h,Nt,Nr,L);          
    %save Y.mat Y
    %save H.mat H
    return
end
%----Calculate MIMO DFE taps based on channel estimation
if(0)
    [C b] = MMSE_MIMO_DFE_LE_Coefficients_fn(h, Nt, Nr, L, K1, K2, K3, SNR);
    N= Kp1;    
    [S1_Equ, Psi] = BLOCK_1_MIMO_MMSE_DFE_fn(Y, QPSK_TxSym, N, C, b, Nt, Nr, K1, K2, K3, Ngsize, Nblk, Nbps);
    S1_Equ= [QPSK_TxSym(:,1:Kp1) S1_Equ];
else
    N= Kp1; %last pilot symbol
    N1= N-K2-L+2; N2= N+Nblk+K1;
    Blk_LLRa_cod= LLRa_cod(:,Nbps*N1-1: Nbps*N2); 
    [x_EstA PsiA LLRe_codA]= MIMO_MMSE_LC_TurboEqu_QPSK_fn(Y,QPSK_TxSym,N,Blk_LLRa_cod,h,L,Nt,Nr,Nblk,Ngsize,Nbps,SNRdB);
end
x_Est= [x_Est x_EstA];
Psi=[Psi PsiA];
LLRe_cod= [LLRe_cod LLRe_codA];
    
%return

BitErr_Num_temp=[];
BER_blk=[];
%------------ for nd-th block equalization ---------------------------
for nd=2:floor((Ndata-Kp1)/Nblk)
    nd
    N = (nd-1)*Nblk+Kp1;   %The last symbol which has been equalized and detected

    %channel re-estimation using detected (or pilot) symbols
    y=Y(:,N-Kp2+1:N); 
    if mod(nd-1,pcyc) == 0
      x= QPSK_TxSym(:,N-Kp2+1:N);
    else
      if Nbps == 1
          x=bpsk(demod_bpsk(x_Est(:,end-Kp2+1:end)));
      elseif Nbps == 2
          x=qpsk(demod_qpsk(x_Est(:,end-Kp2+1:end))); %S1_Equ(N-Np) is synchronized with Y(:,N)
      elseif Nbps == 3
          x=eightpsk(demod_8psk(x_Est(:,end-Kp2+1:end)));    
      elseif Nbps == 4
          x=mqam(demod_mqam(x_Est(:,end-Kp2+1:end),16),16);  
      elseif Nbps == 6
          x=mqam(demod_mqam(x_Est(:,end-Kp2+1:end),64),64);      
      end
    end
    h = TimeDomain_MIMO_ChnnEst_fn(y, x, Nt, Nr, L, Kp2, SNR); %re-estimated channel

      N1= N-K2-L+2; N2= N+Nblk+K1;
      Blk_LLRa_cod= LLRa_cod(:,Nbps*N1-1: Nbps*N2); 
      [x_EstA PsiA LLRe_codA]= MIMO_MMSE_LC_TurboEqu_QPSK_fn(Y,QPSK_TxSym,N,Blk_LLRa_cod,h,L,Nt,Nr,Nblk,Ngsize,Nbps,SNRdB);

    x_Est= [x_Est x_EstA]; %collect new equalized symbols
    Psi=[Psi PsiA+repmat(Psi(:,end),1,length(PsiA))]; %collect accumulated phase estimation
    LLRe_cod= [LLRe_cod LLRe_codA]; %collect LLR for new equalized block

    BitErr_Num_temp=[BitErr_Num_temp sum((demod_qpsk(QPSK_TxSym(:,N+1: N+Nblk))~=demod_qpsk(x_EstA)).').']
    BER_blk= [BER_blk BitErr_Num_temp(:,nd-1)/Nblk/Nbps]
end
    
%------------ for last block equalization ---------------------------
if(mod(Ndata-Kp1,Nblk))
    N = nd*Nblk+Kp1;   %The last symbol which has been equalized and detected
    y=Y(:,N-Kp2+1:N); 

    if mod(nd-1,pcyc) == 0
        x= QPSK_TxSym(:,N-Kp2+1:N);
    else
        if Nbps == 1
            s=bpsk(demod_bpsk(x_Est(:,end-Kp2+1:end)));
        elseif Nbps == 2
            s=qpsk(demod_qpsk(x_Est(:,end-Kp2+1:end))); %S1_Equ(N-Np) is synchronized with Y(:,N)
        elseif Nbps == 3
            s=eightpsk(demod_8psk(x_Est(:,end-Kp2+1:end)));    
        elseif Nbps == 4
            s=mqam(demod_mqam(x_Est(:,end-Kp2+1:end),16),16);  
        elseif Nbps == 6
            s=mqam(demod_mqam(x_Est(:,end-Kp2+1:end),64),64);      
        end
    end
    h = TimeDomain_MIMO_ChnnEst_fn(y, s, Nt, Nr, L, Kp2, SNR); %re-estimated channel

    if(0)
        [C b] = MMSE_MIMO_DFE_LE_Coefficients_fn(h, Nt, Nr, L, K1, K2, K3, SNR);
        [S1_EquA, PsiA] = BLOCK_L_MIMO_MMSE_DFE_fn(Y, QPSK_TxSym, N, C, b, Nt, Nr, K1, K2, K3, Ngsize, Ndata-N, Nbps);
    else
        %Blk_LLRa_cod= LLRa_cod(:,Nbps*N+1: Nbps*(N+Nblk+KK-1)); 
        N1= N-K2-L+2; N2= N+Nblk+K1;
        Blk_LLRa_cod= LLRa_cod(:,Nbps*N1-1: Nbps*N2); 
        [x_EstA PsiA LLRe_codA]= MIMO_MMSE_LC_TurboEqu_QPSK_fn(Y,QPSK_TxSym,N,Blk_LLRa_cod,h,L,Nt,Nr,Nblk,Ngsize,Nbps,SNRdB);
    end

    x_Est= [x_Est x_EstA]; %collect new equalized symbols
    Psi=[Psi PsiA+repmat(Psi(:,end),1,length(PsiA))]; %collect new estimated phase
    LLRe_cod= [LLRe_cod LLRe_codA]; %collect LLR for new equalized block
end

return
