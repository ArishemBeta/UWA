function [doppler_scale,cfo]=D2SearchQPSK(ald,ahd,alc,ahc,P,u,v,sc_idx,Nt,Nr,Ns,K,L,RX_block,pilot_symbol,block_symbol,SNR,SNRdB,Nbps,B)
d=[ald:0.00001:ahd];
c=[alc:0.4:ahc];
% c=alc;(ahc-alc)/2
ld=length(d);
lc=length(c);
recy=zeros(ld,lc);
recy_min=zeros(ld);
nnnn=1;
color=['k.','r.','g.','b.','c.','m.','y.'];
for i=1:ld
    for j=1:lc
        ad=d(i);
        ac=c(j);
%         recx(1,i)=ad;
%         recx(2,j)=ac;
        BER_cost=BER_cost_d_calculation(RX_block,ad,ac,Nt,Nr,block_symbol,Ns,K,L,sc_idx,pilot_symbol,SNR,SNRdB,Nbps,B);
        recy(i,j)=BER_cost;
        scatter3(ad,ac,BER_cost,color(1+mod(i,6)));
        hold on;

        nnnn=nnnn+1;
    end
    recy_min(i)=min(recy(i,:));
end
xlabel('doppler');ylabel('cfo');zlabel('cost');

% image(d,c,recy,'CDataMapping','scaled');%
% colorbar;
% colormap('jet');
% caxis('auto');
% xlabel('doppler','FontSize',12);
% ylabel('cfo(Hz)','FontSize',12);

m=min(min(recy));
% figure();
% plot(d,recy_min);
[id ic]=find(recy==m);
op_doppler_scale=d(id);
op_cfo=c(ic);
doppler_scale=op_doppler_scale;cfo=op_cfo;
return

function BER_cost_d=BER_cost_d_calculation(RX_block,ad,ac,Nt,Nr,block_symbol,Ns,K,L,sc_idx,pilot_symbol,SNR,SNRdB,Nbps,B)
for i=1:Nr
    RX_block(i,:)=interp1((0:length(RX_block(i,:))-1),RX_block(i,:),(0:length(RX_block(i,:))-1)/(1+ad),'spline');
end
for i=1:Nr
    for n=1:length(RX_block(i,:))
        RX_block(i,n)=RX_block(i,n)*(exp(sqrt(-1)*2*pi*(ac-13000*ad)*(n-1)/(Ns*B)));%48828.125
    end
end
y=RX_block(:,1: Ns: K*Ns);
ola=RX_block(:,1+K*Ns: Ns: K*Ns+(L-1)*Ns);
y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
h=FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, Ns, K, L, sc_idx, SNR);
LLa_cod=log(0.5)*ones(2*Nt,Nbps*K);
LLR_info=zeros(Nt,K);
[S_Est,LLe_cod]=MIMO_OFDM_SoftEqu_qpsk_fn_28may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB);
BER_cost_d=0;
symbol_est=S_Est;
pilot_symbol_est=symbol_est(:,sc_idx);
symbol_err=pilot_symbol-pilot_symbol_est;
for nt=1:Nt
    for k=1:length(sc_idx)
        BER_cost_d=BER_cost_d+(symbol_err(nt,k)*conj(symbol_err(nt,k)));
    end
end
return