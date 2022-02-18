 function cfo_iter=CFOIteration(al,ah,P,u,v,sc_idx,Nt,Nr,Ns,K,L,RX_block,pilot_symbol,block_symbol,pilot_bit,block_bit)
% E=0.1;
% maxnum=5;%最大迭代次数
% narvs=1;%目标函数的自变量个数
% particlesize=4;%粒子群规模
% c1=2;%每个粒子的个体学习因子，加速度常数
% c2=2;%每个粒子的社会学习因子，加速度常数
% w=0.6;%惯性因子
% vmax=6;%粒子的最大飞翔速度
% v=2*rand(particlesize,narvs);%粒子飞翔速度
% x=10+20*rand(particlesize,narvs);%粒子所在位置
% %定义适应度函数
% % fitness=inline('(x(1)^2+x(2)^2)/10000','x');
% for i=1:particlesize
% 	f(i)=BER_cost_d_calculation(RX_block,0,x(i),pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol);
% end
% personalbest_x=x;
% personalbest_faval=f;
% [globalbest_faval,i]=min(personalbest_faval);
% globalbest_x=personalbest_x(i,:); 
% k=1;
% while (k<=maxnum)
% 	for i=1:particlesize
% 			f(i)=BER_cost_d_calculation(RX_block,0,x(i,:),pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol);
% 		if f(i)<personalbest_faval(i)
% 			personalbest_faval(i)=f(i);
% 			personalbest_x(i,:)=x(i,:);
% 		end
% 	end
% 	[globalbest_faval,i]=min(personalbest_faval);
% 	globalbest_x=personalbest_x(i,:);
% 	for i=1:particlesize
% 		v(i,:)=w*v(i,:)+c1*rand*(personalbest_x(i,:)-x(i,:))...
% 			+c2*rand*(globalbest_x-x(i,:));
% 		for j=1:narvs
% 			if v(i,j)>vmax
% 				v(i,j)=vmax;
% 			elseif v(i,j)<-vmax
% 				v(i,j)=-vmax;
%             end
% 		end
% 		x(i,:)=x(i,:)+v(i,:);
%     end
%     ff(k)=globalbest_faval;
%     if globalbest_faval<E
%         break
%     end
% %       figure(1)
% %       for i= 1:particlesize
% %       plot(x(i,1),x(i,2),'*')
% %       end
% 	k=k+1;
% end
% cfo_iter=globalbest_x;
% figure(2)
% set(gcf,'color','white');
% plot(1:length(ff),ff)
grad=zeros(1,25);
nnnn=1;
kkkk=0;
delta_l=(ah-al)/5;
delta=(ah-al)/5;
u=1;
P=3;
cfo_cost_rec=zeros(2,25);
BER_cost=0;    
BER_cost=BER_cost_d_calculation(RX_block,0,al,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol);
cfo_cost_rec(1,nnnn)=al;
cfo_cost_rec(2,nnnn)=BER_cost;
while(1)
    if(nnnn>1 && abs(cfo_cost_rec(2,nnnn-1)-cfo_cost_rec(2,nnnn))<0.1)break;end
    a=cfo_cost_rec(1,nnnn);
    BER_cost=cfo_cost_rec(2,nnnn);
    BER_cost_d=BER_cost_d_calculation(RX_block,delta,a,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol);
    grad(nnnn)=2*sqrt(BER_cost)*(sqrt(BER_cost_d)-sqrt(BER_cost))/delta;
    if(nnnn<u) 
        ave_grad=sum(grad(1:nnnn))/nnnn;
    else
        ave_grad=sum(grad(nnnn-u+1:nnnn))/u;
    end
    if(ave_grad<=0)
%         BER_cost=BER_cost+ave_grad*delta;
        if(nnnn>1 && ave_grad>=0.4*grad(nnnn-1) && delta==delta_l && delta~=0.01)%
            kkkk=kkkk+1;
            if(kkkk<P)
                delta_l=delta;
%                 delta=10^(-kkkk);
                delta=delta/10;
            else
                delta=0.01;
            end
            continue;
        end
        BER_cost=BER_cost_d;
        a=a+delta;
        delta_l=delta;
        if(a>ah)
            break;
        end
        nnnn=nnnn+1;
        cfo_cost_rec(1,nnnn)=a;
        cfo_cost_rec(2,nnnn)=BER_cost;
        continue;
    else
        kkkk=kkkk+1;
        if(kkkk<P) 
            delta_l=delta;
%             delta=10^(-kkkk);
            delta=delta/10;
        else
            delta=0.01;
        end
    end
end
% [op_cost op_cost_idx]=min(cfo_cost_rec(2,nnnnn-u+1:nnnnn));
[op_cost op_cost_idx]=min(cfo_cost_rec(2,1:nnnn));
op_cfo=cfo_cost_rec(1,op_cost_idx);
plot(cfo_cost_rec(2,1:nnnn));
hold on;
cfo_iter=op_cfo;
return

function BER_cost_d=BER_cost_d_calculation(RX_block,delta,a,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol)
for i=1:Nr
    for n=1:length(RX_block(i,:))
        RX_block(i,n)=RX_block(i,n)*(exp(sqrt(-1)*2*pi*(a+delta)*(n-1)/48828.125));%39062.5
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