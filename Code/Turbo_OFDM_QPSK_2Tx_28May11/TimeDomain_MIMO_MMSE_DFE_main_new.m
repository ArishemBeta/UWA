%------------------------------------------------
% Time Domain MIMO DFE and LE
% 
%------------------------------------------------

clear all; close all; %pack
rand('state',0); randn('state',0)

addpath ..

%Modulation can be QPSK, 8PSK, 16QAM or 64 QAM
%Nbps =[1 2 3 4 6]; number of bits per symbols
Nt=2;        %number of transmitters
Nr=4;        %number of receivers
Np= 40;      %number of pilots for channel estimation
Ngap= 20;    %number of symbols for zeros (guard time)
Nsym = 9900; %9800;   %number of data symbols
L=11; %11;        %number of channel taps

K2=30;   %Causal part of Feedforward filter (K2 >= 0)
K1=30;  %Anti-causal part of Feedforward filter
K3=0;  %Feedback part, K3 <= K2+L-1
K=K1+K2+1;

Ts=1/4e3;           %symbol interval (seconds)
SNRdB= 30;          %signal-to-noise ratio in dB
SNR= 10^(SNRdB/10); %SNR in linear scale

%--------------------------------------------------------------------
% Load Transmitted and Received Data  
%--------------------------------------------------------------------
if(Nt==1)
    load TimeDomain_SIMO_Tx_Data.mat    %Data_Packet Nbps
    load TimeDomain_SIMO_Rx_Data.mat    %R
    load SIMO_Thetafk.mat               %theta fk
else
    load TimeDomain_MIMO_Tx_Data.mat    %Data_Packet Nbps
    load TimeDomain_MIMO_Rx_Data.mat    %R
    load MIMO_Thetafk.mat               %theta fk
end

N_packet=length(Data_Packet);
Nstart=30+Ngap;
    
Rx=R(:,Nstart+1:N_packet);            %remove the pilots and gap
S1=Data_Packet(:,Nstart+1:N_packet);  %remove the pilots and gap

if(0)
Rx(1,:)=R(1,Nstart+1:N_packet-10);
Rx(2,:)=R(2,Nstart+4:N_packet-7);
Rx(3,:)=R(3,Nstart+3:N_packet-8);
Rx(4,:)=R(4,Nstart+11:N_packet);
end


%------------ for first block equalization ---------------------------
%----------------------------------------------------------
% Channel Estimation for h
%----------------------------------------------------------
N=Np;  %The last symbol of pilots
y=Rx(:,N-Np+1:N); 
s=S1(:,N-Np+1:N);
h = TimeDomain_MIMO_ChnnEst_fn(y, s, Nt, Nr, L, Np, SNR);
if(0)
    load MIMO_Channels.mat
    for m=1: Nr
        figure(100+m)
        clf
        plot(abs(h(m,:)),'-ro'); 
        hold on;
        plot(abs(A(m,:)),'-b*'); 
        legend('Estimated channel','Actual channel');
    end
end

%----------------------------------------------------------
% Channel Equalizer Coefficients, C and b
%----------------------------------------------------------
%K2=30; K1=30; K3=10;  %%K3 <= K2+L-1 
%K=K1+K2+1;
[C b] = MMSE_MIMO_DFE_LE_Coefficients_fn(h, Nt, Nr, L, K1, K2, K3, SNR);

%-----------------------------
% Block channel equalization
%-----------------------------
Ngsize= 10; %10; %4;
Nblock= 40; %500; %500;
if round(Nblock/Ngsize)*Ngsize ~= Nblock
    error('Nblock/Ngsize must be integer')
end
[S1_Equ, Psi] = BLOCK_1_MIMO_MMSE_DFE_fn(Rx, S1, N, C, b, Nt, Nr, K1, K2, K3, Ngsize, Nblock, Nbps);
%S1_Equ: Nt x Nblock;  %Psi: Nt x Nblock;

%------------ for nd-th block equalization ---------------------------
for nd=2:floor((N_packet-Nstart)/Nblock)-2
    nd
  N = (nd-1)*Nblock+Np;   %The last symbol which has been equalized and detected
  y=Rx(:,N-Np+1:N); 

  if Nbps == 1
      s=bpsk(demod_bpsk(S1_Equ(:,N-Np+1-Np:N-Np)));
  elseif Nbps == 2
      s=qpsk(demod_qpsk(S1_Equ(:,N-Np+1-Np:N-Np))); %S1_Equ(N-Np) is synchronized with Rx(:,N)
  elseif Nbps == 3
      s=eightpsk(demod_8psk(S1_Equ(:,N-Np+1-Np:N-Np)));    
  elseif Nbps == 4
      s=mqam(demod_mqam(S1_Equ(:,N-Np+1-Np:N-Np),16),16);  
  elseif Nbps == 6
      s=mqam(demod_mqam(S1_Equ(:,N-Np+1-Np:N-Np),64),64);      
  end

  %s=S1_Equ(1,N-Np+1-Np:N-Np);  %S1_Equ(N-Np) is synchronized with Rx(:,N)
  if mod(nd-1,10) == 0
    s=S1(:,N-Np+1:N);
  end
  h = TimeDomain_MIMO_ChnnEst_fn(y, s, Nt, Nr, L, Np, SNR);

  %figure(100+nd)
  %clf
  %plot(abs(h(1,:)))
  
  [C b] = MMSE_MIMO_DFE_LE_Coefficients_fn(h, Nt, Nr, L, K1, K2, K3, SNR);

  %Ngsize=20;
  %Nblock=500;
  %if mod(nd-1,10) == 0
  %  [S1_EquA, PsiA] = BLOCK_1_MIMO_MMSE_DFE_fn(Rx, S1, N, C, b, Nt, Nr, K1, K2, K3, Ngsize, Nblock, Nbps,Psi(:,length(Psi)));
  %else
    [S1_EquA, PsiA] = BLOCK_n_MIMO_MMSE_DFE_fn(Rx, S1, N, C, b, Nt, Nr, K1, K2, K3, Ngsize, Nblock, Nbps,Psi(:,length(Psi)));    
  %end
  
  S1_Equ=[S1_Equ S1_EquA];
  %Psi=[Psi PsiA+Psi(:,length(Psi))]; %?
  Psi=[Psi PsiA+repmat(Psi(:,end),1,length(PsiA))];
  %X=[X XA];
end

% if K3 == 0
%     if Nr == 1
%         scatterplot(X) %X contains LE-equalized symbols of all receivers
%     else
%         scatterplot(sum(X))
%     end
%     title('Equalized non-Phase-Corrected Symbols')
% end
for n= 1: Nt
    scatterplot(S1_Equ(n,:))
    title('Equalizated and Phase-Corrected Symbols')
end

%S1_Equ(k) is synchronized with S1(Np+k)
S1=S1(:,Np+1:Np+length(S1_Equ(1,:))); 
for n=1: Nt
    if Nbps == 1
        S1_bits(n,:)=demod_bpsk(S1(n,:));
        S1_Equ_bits(n,:)=demod_bpsk(S1_Equ(n,:));     
    elseif Nbps == 2
        S1_bits(n,:)=demod_qpsk(S1(n,:));
        S1_Equ_bits(n,:)=demod_qpsk(S1_Equ(n,:));   
    elseif Nbps == 3
        S1_bits(n,:)=demod_8psk(S1(n,:));
        S1_Equ_bits(n,:)=demod_8psk(S1_Equ(n,:));   
    elseif Nbps == 4
        S1_bits(n,:)=demod_mqam(S1(n,:),16);
        S1_Equ_bits(n,:)=demod_mqam(S1_Equ(n,:),16);   
    elseif Nbps == 6
        S1_bits(n,:)=demod_mqam(S1(n,:),64);
        S1_Equ_bits(n,:)=demod_mqam(S1_Equ(n,:),64);       
    else
        error('Nbps must be 1, 2, 3, 4 or 6 for now')
    end
    Error_Plot(n,:)=[S1_bits(n,:)-S1_Equ_bits(n,:)];    
    BER(n)=sum(abs(S1_bits(n,:)-S1_Equ_bits(n,:)))/length(S1_bits(n,:))
end

%return

for n= 1: Nt    
    figure(11+n)
    clf
    plot(Error_Plot(n,:))
    title('BER')
    xlabel('bit index')
end

for n= 1: Nt  
    figure(200+n)
    clf
    plot(-Psi(n,:))
    title('Phase')
    ylabel('Radian')
    xlabel('group index')
    grid on
end


load MIMO_Thetafk.mat
for n= 1: Nt
    figure(300+n)
    clf
    plot(theta(1,(n-1)*N_packet+1:n*N_packet),'--b','LineWidth',2)
    hold on
    plot(theta(2,(n-1)*N_packet+1:n*N_packet),'-r','LineWidth',2)
    plot(theta(3,(n-1)*N_packet+1:n*N_packet),'--k','LineWidth',1)
    plot(theta(4,(n-1)*N_packet+1:n*N_packet),'-.m','LineWidth',2)
    legend('1','2','3','4')
    ylabel('theta')
    xlabel('symbol index')
    grid on
end

if(1)
    for n= 1: Nt
        figure(310+n)
        clf
        plot(fk(1,(n-1)*N_packet+1:n*N_packet),'--b','LineWidth',2)
        hold on
        plot(fk(2,(n-1)*N_packet+1:n*N_packet),'-r','LineWidth',2)
        plot(fk(3,(n-1)*N_packet+1:n*N_packet),'--k','LineWidth',1)
        plot(fk(4,(n-1)*N_packet+1:n*N_packet),'-.m','LineWidth',2)
        legend('1','2','3','4')
        ylabel('fk')
        xlabel('symbol index')
        grid on
    end
end
