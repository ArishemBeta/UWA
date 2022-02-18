function [doppler_scalee,cfoe,shifte]=D3Iteration(ald,ahd,alc,ahc,P,u,v,sc_idx,Nt,Nr,Ns,K,L,RX_data,Gap,nblk,doppler_scale,pilot_symbol,block_symbol,pilot_bit,block_bit)
grad=zeros(2,100);%-0.0015,-0.0010,18,22
normld=0;normhd=100;normlc=0;normhc=100;
nad=(ahd-ald)/100;nbd=ald;nac=(ahc-alc)/100;nbc=alc;
nnnn=1;
kkkkd=0;kkkkc=0;
delta_lastd=20;delta_lastc=20;%delta_lastd=0.0001 delta_lastc=1
deltad=delta_lastd;deltac=delta_lastc;
ud=1;uc=1;Pd=3;Pc=3;P=1;
rec=zeros(4,100);
BER_cost=0;    
shiftl=-40;shifth=0;
Gapp=Gap;
global gshift;
if(nblk>1)
    for i=2:nblk
        Gapp= Gapp+fix(Ns*1354/(1+doppler_scale(i-1)));%-6770
    end%jump over the pilot block
    RX_block=RX_data(:,Gapp:fix(Gapp+(K*Ns+(L)*Ns)/(1+doppler_scale(nblk-1)))-1);%
else
    RX_block=RX_data(:,Gapp:Gapp+K*Ns+(L)*Ns-1);
end

BER_cost=BER_cost_d_calculation3(RX_block,0,normld,0,normlc,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,nad,nbd,nac,nbc);
rec(1,nnnn)=normld;rec(2,nnnn)=normlc;rec(4,nnnn)=BER_cost;
suc=0;
while(1)
    jumpd=1;jumpc=1;
    normld=rec(1,nnnn);
    normlc=rec(2,nnnn);
    BER_cost=rec(4,nnnn);
    BER_cost_dd=BER_cost_d_calculation3(RX_block,deltad,normld,0,normlc,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,nad,nbd,nac,nbc);
    grad(1,nnnn)=2*sqrt(BER_cost)*(sqrt(BER_cost_dd)-sqrt(BER_cost))/(deltad);%*nad
    BER_cost_dc=BER_cost_d_calculation3(RX_block,0,normld,deltac,normlc,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,nad,nbd,nac,nbc);
    grad(2,nnnn)=2*sqrt(BER_cost)*(sqrt(BER_cost_dc)-sqrt(BER_cost))/(deltac);%*nac
    ave_gradd=sum(grad(1,nnnn));
    ave_gradc=sum(grad(2,nnnn));
    if(ave_gradd<=0)
        if(nnnn>1 && ave_gradd>=0.4*grad(1,nnnn-1) && deltad==delta_lastd && deltad~=10^(-P))%
            kkkkd=kkkkd+1;
            delta_lastd=deltad;
            deltad=max(deltad*0.1,10^(-P));
            jumpd=0;
            suc=0;
        end
    else
        kkkkd=kkkkd+1;
        delta_lastd=deltad;
        deltad=max(deltad*0.1,10^(-P));
        if(kkkkd<Pd) jumpd=0;suc=0;else break;end
    end
    if(ave_gradc<=0)
        if(nnnn>1 && ave_gradc>=0.4*grad(2,nnnn-1) && deltac==delta_lastc && deltac~=10^(-P))%
            kkkkc=kkkkc+1;
            delta_lastc=deltac;
            deltac=max(deltac*0.1,10^(-P));
            jumpc=0;
            suc=0;
        end
    else
        kkkkc=kkkkc+1;
        delta_lastc=deltac;
        deltac=max(deltac*0.1,10^(-P));
        if(kkkkc<Pc) jumpc=0;suc=0;else break;end
    end
%     if(kkkkd>Pd && kkkkc>Pc)break;end
    if((jumpd+jumpc)~=0)
        if(jumpd~=0)normld=normld+deltad;delta_lastd=deltad;end
        if(jumpc~=0)normlc=normlc+deltac;delta_lastc=deltac;end
        BER_cost=BER_cost_d_calculation3(RX_block,0,normld,0,normlc,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,nad,nbd,nac,nbc);
        if(normld>normhd || normlc>normhc)break;end
        nnnn=nnnn+1;suc=suc+1;
        if(suc>=3)deltad=1.2*deltad;deltac=1.2*deltac;else if(suc>=2)deltad=0.1+deltad;deltac=0.1+deltac;end;end
        rec(1,nnnn)=normld;
        rec(2,nnnn)=normlc;
        rec(4,nnnn)=BER_cost;
        scatter3(normld,normlc,BER_cost);
        hold on;
        continue;
    else
        continue;
    end
end
[op_cost op_cost_idx]=min(rec(4,1:nnnn));
doppler_scale(nblk)=rec(1,op_cost_idx)*nad+nbd;
doppler_scalee=doppler_scale(nblk);
cfo(nblk)=rec(2,op_cost_idx)*nac+nbc;
cfoe=cfo(nblk);
shifte=0;
% plot(rec(4,1:nnnn),'.');
gshift=-40;
% while(1)
%     Gapp=Gapp+gshift;
%     if(nblk>1)
%         RX_block=RX_data(:,Gapp:fix((Gapp+K*Ns+(L)*Ns)/(1+doppler_scale(nblk-1)))-1);
%     else
%         RX_block=RX_data(:,Gapp:Gapp+K*Ns+(L)*Ns-1);
%     end
%     for i=1:Nr
%         RX_block(i,:)=interp1((0:length(RX_block(i,:))-1),RX_block(i,:),(0:length(RX_block(i,:))-1)/(1+doppler_scale(nblk)),'spline');
%     end
%     for i=1:Nr
%         for n=1:length(RX_block(i,:))
%             RX_block(i,n)=RX_block(i,n)*(exp(sqrt(-1)*2*pi*(cfo(nblk)*(n-1)/48828.125)));%39062.5
%         end
%     end
%     y=RX_block(:,1: Ns: K*Ns);
%     ola=RX_block(:,1+K*Ns: Ns: K*Ns+(L-1)*Ns);
%     y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
%     h=FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, Ns, K, L, sc_idx, 10);
%     LLa_cod= log(0.5)*ones(2*Nt,4096);
%     LLR_info= zeros(Nt,2048);
%     [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_16qam_fn_31may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,4,10);
%     symbol_est=S_Est;
%     pilot_symbol_est=symbol_est(:,sc_idx);
%     symbol_err=pilot_symbol-pilot_symbol_est;
%     BER_cost_d=0;
%     for nt=1:Nt
%         for k=1:K/4
%             BER_cost_d=BER_cost_d+symbol_err(nt,k)*conj(symbol_err(nt,k));
%         end
%     end
% %     if(BER_cost_d<=op_cost)
% %         op_cost=BER_cost_d;
% %         gshift=gshift+5;
% %     else
% %         gshift=gshift-5;
% %         shifte=gshift;
% %         break;
% %     end
%     gshift=gshift+5;
%     plot(gshift,BER_cost_d,'*');hold on;
%     if(gshift>40)break;end
%     
%     end
return

function BER_cost_d=BER_cost_d_calculation3(RX_block,deltad,ad,deltac,ac,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,nad,nbd,nac,nbc)
for i=1:Nr
    RX_block(i,:)=interp1((0:length(RX_block(i,:))-1),RX_block(i,:),(0:length(RX_block(i,:))-1)/(1+(ad+deltad)*nad+nbd),'spline');
end
for i=1:Nr
    for n=1:length(RX_block(i,:))
        RX_block(i,n)=RX_block(i,n)*(exp(sqrt(-1)*2*pi*((ac+deltac)*nac+nbc)*(n-1)/48828.125));%39062.5
    end
end
y=RX_block(:,1: Ns: K*Ns);
ola=RX_block(:,1+K*Ns: Ns: K*Ns+(L-1)*Ns);
y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
h=FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, Ns, K, L, sc_idx, 10);
LLa_cod= log(0.5)*ones(2*Nt,4096);
LLR_info= zeros(Nt,2048);
[S_Est LLe_cod]= MIMO_OFDM_SoftEqu_16qam_fn_31may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,4,10);
BER_cost_d=0;
symbol_est=S_Est;
pilot_symbol_est=symbol_est(:,sc_idx);
symbol_err=pilot_symbol-pilot_symbol_est;
for nt=1:Nt
    for k=1:K/4
        BER_cost_d=BER_cost_d+symbol_err(nt,k)*conj(symbol_err(nt,k));
    end
end
return