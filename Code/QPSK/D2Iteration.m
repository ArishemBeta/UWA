function [doppler_scale,cfo]=D2Iteration(ald,ahd,alc,ahc,P,u,v,sc_idx,Nt,Nr,Ns,K,L,RX_block,pilot_symbol,block_symbol,pilot_bit,block_bit)
grad=zeros(2,100);%-0.0015,-0.0010,18,22
normld=0;normhd=100;normlc=0;normhc=100;
nad=(ahd-ald)/100;nbd=ald;nac=(ahc-alc)/100;nbc=alc;
nnnn=1;
kkkkd=0;kkkkc=0;
% delta_lastd=(ahd-ald)/5;delta_lastc=round((ahc-alc)/5);%delta_lastd=0.0001 delta_lastc=1
% deltad=delta_lastd;deltac=delta_lastc;
delta_lastd=20;delta_lastc=20;%delta_lastd=0.0001 delta_lastc=1
deltad=delta_lastd;deltac=delta_lastc;
ud=1;uc=1;Pd=3;Pc=3;P=1;
rec=zeros(4,100);
BER_cost=0;    
BER_cost=BER_cost_d_calculation(RX_block,0,normld,0,normlc,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,nad,nbd,nac,nbc);
rec(1,nnnn)=normld;rec(2,nnnn)=normlc;rec(4,nnnn)=BER_cost;
suc=0;
while(1)
    jumpd=1;jumpc=1;
    normld=rec(1,nnnn);
    normlc=rec(2,nnnn);
    BER_cost=rec(4,nnnn);
    BER_cost_dd=BER_cost_d_calculation(RX_block,deltad,normld,0,normlc,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,nad,nbd,nac,nbc);
    grad(1,nnnn)=2*sqrt(BER_cost)*(sqrt(BER_cost_dd)-sqrt(BER_cost))/(deltad);%*nad
    BER_cost_dc=BER_cost_d_calculation(RX_block,0,normld,deltac,normlc,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,nad,nbd,nac,nbc);
    grad(2,nnnn)=2*sqrt(BER_cost)*(sqrt(BER_cost_dc)-sqrt(BER_cost))/(deltac);%*nac
    ave_gradd=sum(grad(1,nnnn));
    ave_gradc=sum(grad(2,nnnn));
    if(ave_gradd<=0)
%         if(nnnn>1 && ave_gradd>=0.9*grad(1,nnnn-1) && deltad==delta_lastd && deltad~=10^(-P))%
%             kkkkd=kkkkd+1;
%             delta_lastd=deltad;
%             deltad=max(deltad*0.5,10^(-P));
%             jumpd=0;
%             suc=0;
%         end
    else
        kkkkd=kkkkd+1;
        delta_lastd=deltad;
        deltad=max(deltad*0.5,10^(-P));
        if(kkkkd<Pd) jumpd=0;suc=0;else break;end
    end
    if(ave_gradc<=0)
%         if(nnnn>1 && ave_gradc>=0.9*grad(2,nnnn-1) && deltac==delta_lastc && deltac~=10^(-P))%
%             kkkkc=kkkkc+1;
%             delta_lastc=deltac;
%             deltac=max(deltac*0.5,10^(-P));
%             jumpc=0;
%             suc=0;
%         end
    else
        kkkkc=kkkkc+1;
        delta_lastc=deltac;
        deltac=max(deltac*0.5,10^(-P));
        if(kkkkc<Pc) jumpc=0;suc=0;else break;end
    end
%     if(kkkkd>Pd && kkkkc>Pc)break;end
    if((jumpd+jumpc)~=0)
        if(jumpd~=0)normld=normld+deltad;delta_lastd=deltad;end
        if(jumpc~=0)normlc=normlc+deltac;delta_lastc=deltac;end
        BER_cost=BER_cost_d_calculation(RX_block,0,normld,0,normlc,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,nad,nbd,nac,nbc);
        if(normld>normhd || normlc>normhc)break;end
        nnnn=nnnn+1;suc=suc+1;
        if(suc>=3)deltad=1.2*deltad;deltac=1.2*deltac;else if(suc>=2)deltad=0.1+deltad;deltac=0.1+deltac;end;end
        rec(1,nnnn)=normld;
        rec(2,nnnn)=normlc;
        rec(4,nnnn)=BER_cost;
%         scatter3(normld,normlc,BER_cost);
%         hold on;
        continue;
    else
        continue;
    end
end
[op_cost op_cost_idx]=min(rec(4,1:nnnn));
doppler_scale=rec(1,op_cost_idx)*nad+nbd;
cfo=rec(2,op_cost_idx)*nac+nbc;
plot(rec(4,1:nnnn),'.');
% grad=zeros(2,100);%-0.0015,-0.0010,18,22
% norml=0;normh=100;
% aa1=(ahd-ald)/100;bb1=ald;
% aa2=(ahc-alc)/100;bb2=alc;
% nnnn=1;
% kkkkd=4;kkkkc=0;
% delta_lastd=(ahd-ald)/5;delta_lastc=round((ahc-alc)/5);%delta_lastd=0.0001 delta_lastc=1
% deltad=delta_lastd;deltac=delta_lastc;
% ud=1;uc=1;
% Pd=6;Pc=2;
% rec=zeros(4,100);
% BER_cost=0;    
% if(nnnn==1)
%     BER_cost=BER_cost_d_calculation(RX_block,0,ald,0,alc,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol);
%     rec(1,nnnn)=ald;
%     rec(2,nnnn)=alc;
%     rec(4,nnnn)=BER_cost;
% end
% while(1)
%     jumpd=1;jumpc=1;
%     ad=rec(1,nnnn);
%     ac=rec(2,nnnn);
%     BER_cost=rec(4,nnnn);
%     BER_cost_dd=BER_cost_d_calculation(RX_block,deltad,ad,0,ac,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol);
%     grad(1,nnnn)=2*sqrt(BER_cost)*(sqrt(BER_cost_dd)-sqrt(BER_cost))/deltad;
%     BER_cost_dc=BER_cost_d_calculation(RX_block,0,ad,deltac,ac,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol);
%     grad(2,nnnn)=2*sqrt(BER_cost)*(sqrt(BER_cost_dc)-sqrt(BER_cost))/deltac;
%     if(nnnn<u) 
%         ave_grad=sum(grad(1:nnnn))/nnnn;
%     else
%         ave_gradd=sum(grad(1,nnnn-u+1:nnnn))/u;
%         ave_gradc=sum(grad(2,nnnn-u+1:nnnn))/u;
%     end
%     if(ave_gradd<=0)
%         if(nnnn>1 && ave_gradd>=0.4*grad(1,nnnn-1) && deltad==delta_lastd && deltad~=10^(-Pd))%
%             kkkkd=kkkkd+1;
%             if(kkkkd<P)
%                 delta_lastd=deltad;
%                 deltad=deltad*0.2;
%             else
%                 deltad=10^(-Pd);
%             end
%             jumpd=0;
% %             continue;
%         else
% %             ad=ad+deltad;
% %             delta_lastd=deltad;
%         end
%     else
%         kkkkd=kkkkd+1;
%         if(kkkkd<Pd) 
%             delta_lastd=deltad;
%             deltad=deltad*0.2;
%         else
%             deltad=10^(-Pd);
%         end
%         if(kkkkd<=Pd) 
%             jumpd=0;
% %             continue;
%         else
%             break;
%         end
%     end
%     if(ave_gradc<=0)
%         if(nnnn>1 && ave_gradc>=0.4*grad(2,nnnn-1) && deltac==delta_lastc && deltac~=10^(-Pc))%
%             kkkkc=kkkkc+1;
%             if(kkkkc<P)
%                 delta_lastc=deltac;
% %                 deltac=10^(-kkkkc);
%                 deltac=deltac*0.2;
%             else
%                 deltac=10^(-Pc);
%             end
%             jumpc=0;
% %             continue;
%         end
%     else
%         kkkkc=kkkkc+1;
%         if(kkkkc<Pc) 
%             delta_lastc=deltac;
% %             deltac=10^(-kkkkc);
%             deltac=deltac*0.2;
%         else
%             deltac=10^(-Pc);
%         end
%         if(kkkkc<=Pc) 
%             jumpc=0;
% %             continue;
%         else
%             break;
%         end
%     end
% %     if(kkkkd>Pd && kkkkc>Pc)break;end
%     if((jumpd+jumpc)~=0)
%         if(jumpd~=0)ad=ad+deltad;delta_lastd=deltad;end
%         if(jumpc~=0)ac=ac+deltac;delta_lastc=deltac;end
%         BER_cost=BER_cost_d_calculation(RX_block,0,ad,0,ac,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol);
%         if(ad>ahd || ac>ahc)break;end
%         nnnn=nnnn+1;
%         rec(1,nnnn)=ad;
%         rec(2,nnnn)=ac;
%         rec(4,nnnn)=BER_cost;
% %         scatter3(ad,ac,BER_cost);
% %         hold on;
%         continue;
%     else
%         continue;
%     end
% end
% % [op_cost op_cost_idx]=min(cfo_cost_rec(2,nnnnn-u+1:nnnnn));
% [op_cost op_cost_idx]=min(rec(4,1:nnnn));
% op_doppler_scale=rec(1,op_cost_idx);
% op_cfo=rec(2,op_cost_idx);
% plot(rec(4,1:nnnn),'.');
% doppler_scale=op_doppler_scale;cfo=op_cfo;
return


% function BER_cost_d=BER_cost_d_calculation(RX_block,deltad,ad,deltac,ac,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol)
% for i=1:Nr
%     RX_block(i,:)=interp1((0:length(RX_block(i,:))-1),RX_block(i,:),(0:length(RX_block(i,:))-1)/(1+ad+deltad),'spline');
% end
% for i=1:Nr
%     for n=1:length(RX_block(i,:))
%         RX_block(i,n)=RX_block(i,n)*(exp(sqrt(-1)*2*pi*(ac+deltac)*(n-1)/48828.125));%39062.5
%     end
% end
% y=RX_block(:,1: Ns: K*Ns);
% ola=RX_block(:,1+K*Ns: Ns: K*Ns+(L-1)*Ns);
% y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
% h=FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, Ns, K, L, sc_idx, 10);
% LLa_cod=log(0.5)*ones(2*Nt,2048);
% LLR_info=zeros(Nt,1024);
% [S_Est,LLe_cod]=MIMO_OFDM_SoftEqu_qpsk_fn_28may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,2,10);
% BER_cost_d=0;
% symbol_est=S_Est;
% pilot_symbol_est=symbol_est(:,sc_idx);
% symbol_err=pilot_symbol-pilot_symbol_est;
% for nt=1:Nt
%     for k=1:K/4
%         BER_cost_d=BER_cost_d+symbol_err(nt,k)*conj(symbol_err(nt,k));
%     end
% end
% return


function BER_cost_d=BER_cost_d_calculation(RX_block,deltad,ad,deltac,ac,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,nad,nbd,nac,nbc)
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
LLa_cod=log(0.5)*ones(2*Nt,2048);
LLR_info=zeros(Nt,1024);
[S_Est,LLe_cod]=MIMO_OFDM_SoftEqu_qpsk_fn_28may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,2,10);
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