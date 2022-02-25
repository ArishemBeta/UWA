function [doppler_scale,cfo]=D2Search(ald,ahd,alc,ahc,P,u,v,sc_idx,Nt,Nr,Ns,K,L,RX_block,pilot_symbol,block_symbol,SNR,SNRdB)
recx=zeros(2,300);
d=[ald:0.0001:ahd];
c=[alc:0.1:ahc];
% c=alc;(ahc-alc)/2
ld=length(d);
lc=length(c);
recy=zeros(ld,lc);
nnnn=1;
color=['k.','r.','g.','b.','c.','m.','y.'];
for i=1:ld
    for j=1:lc
        ad=d(i);
        ac=c(j);
%         recx(1,i)=ad;
%         recx(2,j)=ac;
        BER_cost=BER_cost_d_calculation(RX_block,0,ad,0,ac,Nt,Nr,block_symbol,Ns,K,L,sc_idx,pilot_symbol,SNR,SNRdB);
        recy(i,j)=BER_cost;
        scatter3(ad,ac,BER_cost,color(1+mod(i,6)));
        hold on;
        nnnn=nnnn+1;
    end
end
% plot(d,recy.');
% hold on;
xlabel('doppler');ylabel('cfo');zlabel('cost');
% [dd,cc]=meshgrid(d,c);
% scatter3(d,c,recy);
m=min(min(recy));
[id ic]=find(recy==m);
op_doppler_scale=d(id);
op_cfo=c(ic);
doppler_scale=op_doppler_scale;cfo=op_cfo;
return

function BER_cost_d=BER_cost_d_calculation(RX_block,deltad,ad,deltac,ac,Nt,Nr,block_symbol,Ns,K,L,sc_idx,pilot_symbol,SNR,SNRdB)
for i=1:Nr
    RX_block(i,:)=interp1((0:length(RX_block(i,:))-1),RX_block(i,:),(0:length(RX_block(i,:))-1)/(1+ad+deltad),'spline');
end
for i=1:Nr
    for n=1:length(RX_block(i,:))
        RX_block(i,n)=RX_block(i,n)*(exp(sqrt(-1)*2*pi*(ac+deltac-13000*(ad+deltad))*(n-1)/(Ns*9765.625)));%48828.125
    end
end
y=RX_block(:,1: Ns: K*Ns);
ola=RX_block(:,1+K*Ns: Ns: K*Ns+(L-1)*Ns);
y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
h=FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, Ns, K, L, sc_idx, SNR);
LLa_cod=log(0.5)*ones(2*Nt,2048);
LLR_info=zeros(Nt,1024);
[S_Est,LLe_cod]=MIMO_OFDM_SoftEqu_qpsk_fn_28may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,2,SNRdB);
BER_cost_d=0;
symbol_est=S_Est;
pilot_symbol_est=symbol_est(:,sc_idx);
symbol_err=pilot_symbol-pilot_symbol_est;
for nt=1:Nt
    for k=1:K/4
        BER_cost_d=BER_cost_d+(symbol_err(nt,k)*conj(symbol_err(nt,k)));
    end
end
return