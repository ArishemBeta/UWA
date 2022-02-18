 function cfo_iter=CFOIteration(al,ah,P,u,v,sc_idx,Nt,Nr,Ns,K,L,RX_block,pilot_symbol,block_symbol,pilot_bit,block_bit)
grad=zeros(1,25);
nnnn=1;kkkk=0;
norml=0;normh=100;
delta_l=20;delta=20;
na=(ah-al)/100;nb=al;
u=1;P=3;
cfo_cost_rec=zeros(2,25);
BER_cost=BER_cost_d_calculation(RX_block,0,norml,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,na,nb);
cfo_cost_rec(1,nnnn)=norml;
cfo_cost_rec(2,nnnn)=BER_cost;
suc=0;
while(1)
%     if(nnnn>1 && abs(cfo_cost_rec(2,nnnn-1)-cfo_cost_rec(2,nnnn))<0.1)break;end
    norml=cfo_cost_rec(1,nnnn);
    BER_cost=cfo_cost_rec(2,nnnn);
    BER_cost_d=BER_cost_d_calculation(RX_block,delta,norml,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,na,nb);
    grad(nnnn)=2*sqrt(BER_cost)*(sqrt(BER_cost_d)-sqrt(BER_cost))/delta;
    if(nnnn<u) 
        ave_grad=sum(grad(1:nnnn))/nnnn;
    else
        ave_grad=sum(grad(nnnn-u+1:nnnn))/u;
    end
    if(ave_grad<=0)
%         BER_cost=BER_cost+ave_grad*delta;
        if(nnnn>1 && ave_grad>=0.4*grad(nnnn-1) && delta==delta_l && delta~=10^(-P))%
            kkkk=kkkk+1;
            delta_l=delta;
            delta=max(delta*0.1,10^(-P));
            if(kkkk>=P)break;end
            continue;
        end
        norml=norml+delta;
        delta_l=delta;
        BER_cost=BER_cost_d_calculation(RX_block,0,norml,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,na,nb);
        if(norml>normh)break;end
        nnnn=nnnn+1;
        suc=suc+1;
        if(suc>=3)delta=1.2*delta;else if(suc>=2)delta=0.1+delta;end;end
        cfo_cost_rec(1,nnnn)=norml;
        cfo_cost_rec(2,nnnn)=BER_cost;
        continue;
    else
        kkkk=kkkk+1;
        delta_l=delta;
        delta=max(delta*0.1,10^(-P));
        if(kkkk>=P)break;end
    end
end
% [op_cost op_cost_idx]=min(cfo_cost_rec(2,nnnnn-u+1:nnnnn));
[op_cost op_cost_idx]=min(cfo_cost_rec(2,1:nnnn));
op_cfo=cfo_cost_rec(1,op_cost_idx)*na+nb;
plot(cfo_cost_rec(2,2:nnnn),'*');
cfo_iter=op_cfo;
return

function BER_cost_d=BER_cost_d_calculation(RX_block,delta,a,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol,na,nb)
for i=1:Nr
    for n=1:length(RX_block(i,:))
        RX_block(i,n)=RX_block(i,n)*(exp(sqrt(-1)*2*pi*((a+delta)*na+nb)*(n-1)/48828.125));%39062.5
    end
end
y=RX_block(:,1: Ns: K*Ns);
ola=RX_block(:,1+K*Ns: Ns: K*Ns+(L-1)*Ns);
y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
h=FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, Ns, K, L, sc_idx, 10);
LLa_cod=log(0.5)*ones(2*Nt,4096);
LLR_info=zeros(Nt,2048);
[S_Est,LLe_cod]=MIMO_OFDM_SoftEqu_16qam_fn_31may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,4,10);
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