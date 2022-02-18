%Author: Jun Tao
%Place: MST
%Date: March 02, 2010
%Purpose: MIMO turbo LE for the whole input block
%
%
function [x_Est LLe_cod x_Equ_Soft]= Soft_MIMO_MMSE_LE_packet(Y,X,LLa_cod,h,Ndata,K1,K2,L,Nt,Nr,Nbps,SNRdB);
%function [x_Est LLe_cod]= Soft_MIMO_MMSE_LE_packet(Yblk,QPSK_TxSym_blk,LLa_cod,h,K,K1,K2,L,Nt,Nr,Nbps,SNRdB);

SNR= 10^(SNRdB/10);  %SNR in linear scale
Kp= 400; %pilot length
LLe_cod= []; 
x_Est= [];
x_Equ_Soft= [];
pcyc= 30;

Nsblk= 1024; %equalize the whole block at once
Nrepeat= Ndata/Nsblk; %total number of sub-blocks to be processed
KK= K1+K2+L;

head= 1;
tail= 1;

%----perform soft-decision equalization
LLa_cod_blk= LLa_cod(:,1:(Nsblk+KK-1)*Nbps); 
[x_EstA LLe_codA x_Equ_SoftA]= MIMO_MMSE_LC_TurboEqu_QPSK_fn(Y(:,1:Nsblk+K1+K2),LLa_cod_blk,h,L,Nt,Nr,Nsblk,Nbps,SNRdB,head,tail);


x_Est= [x_Est x_EstA];
x_Equ_Soft= [x_Equ_Soft x_Equ_SoftA];
LLe_cod= [LLe_cod LLe_codA];

BitErr_Num_temp=[];
BER_blk=[];


%------------ for nd-th block equalization ---------------------------
for nd=2: Nrepeat
  nd
  N = (nd-1)*Nsblk;   %The last symbol which has been equalized and detected

  %channel re-estimation using detected (or pilot) symbols
  y=Y(:,N-Kp+1:N); 
  if mod(nd-1,pcyc) == 0
      x= X(:,N-Kp+1:N);
  else
      if Nbps == 1
          x=bpsk(demod_bpsk(x_Est(:,end-Kp+1:end)));
      elseif Nbps == 2
          x=qpsk(demod_qpsk(x_Est(:,end-Kp+1:end)));
      elseif Nbps == 3
          x=eightpsk(demod_8psk(x_Est(:,end-Kp+1:end)));    
      elseif Nbps == 4
          x=mqam(demod_mqam(x_Est(:,end-Kp+1:end),16),16);  
      elseif Nbps == 6
          x=mqam(demod_mqam(x_Est(:,end-Kp+1:end),64),64);      
      end
  end
  h = TimeDomain_MIMO_ChnnEst_fn(y, x, Nt, Nr, L, Kp, SNR); %re-estimated channel

  N1= N-K2-L+2; N2= N+Nsblk+K1;
  LLa_cod_blk= LLa_cod(:,(N1-1)*Nbps+1: Nbps*N2); 
  [x_EstA LLe_codA]= MIMO_MMSE_LC_TurboEqu_QPSK_fn(Y(:,N1+L-1:N2),LLa_cod_blk,h,L,Nt,Nr,Nsblk,Nbps,SNRdB);

  x_Est= [x_Est x_EstA]; %collect new equalized symbols
  LLe_cod= [LLe_cod LLe_codA]; %collect LLR for new equalized block
  if(1)
  LLa_cod(:,Nbps*N+1: Nbps*(N+Nblk))= LLe_codA+LLa_cod(:,Nbps*N+1: Nbps*(N+Nblk)); 
  end
  
  BitErr_Num_temp=[BitErr_Num_temp sum((demod_qpsk(X(:,N+1: N+Nsblk))~=demod_qpsk(x_EstA)).').']
  BER_blk= [BER_blk BitErr_Num_temp(:,nd-1)/Nsblk/Nbps]
end

if(0)
    %------------ for last block equalization ---------------------------
    if(mod(Ndata,Nsblk))
        N = nd*Nsblk+Kp1;   %The last symbol which has been equalized and detected
        y=Y(:,N-Kp+1:N); 

        if mod(nd-1,pcyc) == 0
            x= X(:,N-Kp+1:N);
        else
            if Nbps == 1
                s=bpsk(demod_bpsk(x_Est(:,end-Kp+1:end)));
            elseif Nbps == 2
                s=qpsk(demod_qpsk(x_Est(:,end-Kp+1:end))); %S1_Equ(N-Np) is synchronized with Y(:,N)
            elseif Nbps == 3
                s=eightpsk(demod_8psk(x_Est(:,end-Kp+1:end)));    
            elseif Nbps == 4
                s=mqam(demod_mqam(x_Est(:,end-Kp+1:end),16),16);  
            elseif Nbps == 6
                s=mqam(demod_mqam(x_Est(:,end-K2+1:end),64),64);      
            end
        end
        h = TimeDomain_MIMO_ChnnEst_fn(y, s, Nt, N_chnn, L, Kp, SNR); %re-estimated channel


        N1= N-K2-L+2; N2= N+Nblk+K1;
        LLa_cod_blk= LLa_cod(:,Nbps*N1-(Nbps-1): Nbps*N2); 
        [x_EstA LLe_codA]= MIMO_MMSE_LC_TurboEqu_QPSK_fn(Y(:,N1+L-1:end),LLa_cod_blk,h,L,Nt,Nr,Nblk,Nbps,SNRdB);

        x_Est= [x_Est x_EstA]; %collect new equalized symbols
        LLe_cod= [LLe_cod LLe_codA]; %collect LLR for new equalized block
        if(1)
        LLa_cod(:,Nbps*N+1: Nbps*(N+Nblk))= LLe_codA+LLa_cod(:,Nbps*N+1: Nbps*(N+Nblk));  
        end
    end
end

return
