function [doppler_scalee,cfoe,shifte]=D3Fit(ald,ahd,alc,ahc,sc_idx,Nt,Nr,Ns,K,L,RX_data,nblk,doppler_scale,block_symbol,SNRdB)
normld=0;normhd=100;normlc=0;normhc=100;
nad=(ahd-ald)/100;nbd=ald;nac=(ahc-alc)/100;nbc=alc;
BER_cost=0;    
shiftl=-40;shifth=0;
global gshift;
Kg=240;
Npoint=10;
Noffset=(nblk-1)*(K+Kg);
if(nblk>1)
    for i=2:nblk
        Noffset=Noffset-round(doppler_scale(i-1)*(K+Kg));%start of the nblk-th block
    end
    RX_block=RX_data(:,Noffset+1:Noffset+round((K+Kg)*(1-doppler_scale(nblk-1))));
else
    RX_block=RX_data(:,Noffset+1:Noffset+(K+Kg));
end
for i=1:Npoint
    sampd(i)=(i-1)*100/Npoint+100/Npoint*0.5; %rand()
    BER_cost(i)=cost_fit_cal(RX_block,sampd(i),50,Nt,Nr,block_symbol,Ns,K,L,sc_idx,nad,nbd,nac,nbc,SNRdB);
end
plot(sampd*nad+nbd,BER_cost,'*');
e=polyfit(sampd,BER_cost,4);
Y=polyval(e,[0:1:100]);
hold on;
plot([0:1:100]*nad+nbd,Y,'.');
[a,b]=min(Y);
doppler_scalee=b*nad+nbd;
% cfo(nblk)=rec(2,op_cost_idx)*nac+nbc;
cfoe=0;
shifte=0;
return

function BER_cost_d=cost_fit_cal(RX_block,ad,ac,Nt,Nr,block_symbol,Ns,K,L,sc_idx,nad,nbd,nac,nbc,SNRdB)
[I,D]=numden(sym(rats((ad)*nad+nbd)));
I=eval(I);D=eval(D);
clear Rx_block;
for i=1:Nr
    Rx_block(i,:)=resample(RX_block(i,:),D+I,D);
%     RX_block(i,:)=interp1((0:length(RX_block(i,:))-1),RX_block(i,:),(0:length(RX_block(i,:))-1)/(1+(ad)*nad+nbd),'spline');
end
% for i=1:Nr
%     for n=1:length(RX_block(i,:))
%         RX_block(i,n)=RX_block(i,n)*(exp(sqrt(-1)*2*pi*((ac)*nac+nbc)*(n-1)/9765.625));%39062.5
%     end
% end
y=Rx_block(:,1: Ns: K*Ns);
ola=Rx_block(:,1+K*Ns: Ns: K*Ns+(L-1)*Ns);
y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
h=FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, Ns, K, L, sc_idx, SNRdB);
if(0)
    LLa_cod= log(0.5)*ones(2*Nt,4096);
    LLR_info= zeros(Nt,2048);
    [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_16qam_fn_31may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,4,SNRdB);
else
    LLa_cod=log(0.5)*ones(2*Nt,2048);
    LLR_info=zeros(Nt,1024);
    [S_Est,LLe_cod]=MIMO_OFDM_SoftEqu_qpsk_fn_28may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,2,SNRdB);
end
BER_cost_d=0;
symbol_est=S_Est;
p_idx=[1:4:K];%K*3/8+1:K*5/8 1:K/8,K*7/8+1:K
pilot_symbol=block_symbol(p_idx);
pilot_symbol_est=symbol_est(:,p_idx);
symbol_err=pilot_symbol-pilot_symbol_est;
for nt=1:Nt
    for k=1:K/4
        BER_cost_d=BER_cost_d+(symbol_err(nt,k)*conj(symbol_err(nt,k)));
    end
end
return