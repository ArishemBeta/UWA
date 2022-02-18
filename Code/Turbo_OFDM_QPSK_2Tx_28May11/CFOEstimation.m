function cfo = CFOEstimation(al,ah,P,u,v,sc_idx,Nt,Nr,Ns,K,L,RX_block,pilot_symbol,block_symbol,pilot_bit,block_bit)
% al估计范围下限 ah估计范围上限 v最小分辨率 单位Hz
% P计数阈值 u平均窗口长度 sc_idx导频序号 Nt发射天线 Nr接收天线 Ns过采样因子 K载波数量 QPSK_TxSym发送QPSK符号
% nblk当前块序号 L信道长度 RX_data接收信号 Gap块起始点 Chnn_idx接收机序号
k=0;%计数器
al=-0.06;%-0.06,0.03
ah=0.03;
% a=zeros(Nt,Nr);
% a(:,:)=al;
delta=0.001;
candidate_a=(al:delta:ah);
Niter=length(candidate_a);
candidate_cost=zeros(1,Niter);
for niter=1:Niter
    a=candidate_a(niter);
    BER_cost = cBER_cost_function(RX_block,delta,a,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol);
    candidate_cost(niter)=BER_cost;
end
[min_cost,index]=min(candidate_cost);
cfo=candidate_a(index);
return

function cBER_cost = cBER_cost_function(RX_block,delta,a,pilot_bit,Nt,Nr,block_symbol,Ns,K,L,sc_idx,block_bit,pilot_symbol)
for i=1:Nr
    for n=1:length(RX_block(i,:))
        RX_block(i,n)=RX_block(i,n)*(exp(i*2*pi*a*(n-1)/39062.5));
    end
end
% y=RX_block(:,1441: Ns: K*Ns+1440);
% ola=RX_block(:,1441+K*Ns: Ns: K*Ns+(L-1)*Ns+1440);
y=RX_block(:,1: Ns: K*Ns);
ola=RX_block(:,1+K*Ns: Ns: K*Ns+(L-1)*Ns);
y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
h=FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, Ns, K, L, sc_idx, 10);
LLa_cod=log(0.5)*ones(2*Nt,2048);
LLR_info=zeros(Nt,1024);
[S_Est,LLe_cod]=MIMO_OFDM_SoftEqu_qpsk_fn_28may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,2,10);
cBER_cost=0;
if(1)%Symbol Estimation
    symbol_est=S_Est;
    pilot_symbol_est=symbol_est(:,sc_idx);
%     phase_p(1,:)=angle(pilot_symbol(1,:));
%     phase_p(2,:)=angle(pilot_symbol(2,:));
%     phase_pe(1,:)=angle(pilot_symbol_est(1,:));
%     phase_pe(2,:)=angle(pilot_symbol_est(2,:));
    symbol_err=pilot_symbol-pilot_symbol_est;
%     symbol_err=phase_p-phase_pe;
%     scatterplot(pilot_symbol_est(1,:)); title('Tx 1');
%     scatterplot(pilot_symbol_est(2,:)); title('Tx 2');
    for nt=1:Nt
        for k=1:K/4
            cBER_cost=cBER_cost+symbol_err(nt,k)*conj(symbol_err(nt,k));
        end
    end
else %LLR Estimation
    for nt= 1: Nt
        LLe_cod((nt-1)*2+1,:)= randdeintrlv(LLe_cod((nt-1)*2+1,:), 0);
        LLe_cod((nt-1)*2+2,:)= randdeintrlv(LLe_cod((nt-1)*2+2,:), 0);
        [LLR_info(nt,:),LLe_cod((nt-1)*2+1: nt*2,:)] = MAPConvDecoder_R1(LLe_cod((nt-1)*2+1: nt*2,:),1024);
        LLe_cod((nt-1)*2+1,:)= randintrlv(LLe_cod((nt-1)*2+1,:), 0);
        LLe_cod((nt-1)*2+2,:)= randintrlv(LLe_cod((nt-1)*2+2,:), 0);
        LLa_cod((nt-1)*2+1: nt*2,:)= LLe_cod((nt-1)*2+1: nt*2,:);
    end
    Decode_Bit=zeros(Nt,K);
    Decode_PilotBit=zeros(Nt,K/4);
    BER=zeros(Nt,K/4);
    for nt=1:Nt
        Decode_Bit(nt,:)=LLR_info(nt,:)<0;
        Decode_PilotBit(nt,:)=Decode_Bit(nt,sc_idx);
        BER(nt,:)=xor(pilot_bit(nt,:),Decode_PilotBit(nt,:));
%         BER(nt,:)=xor(block_bit(nt,:),Decode_Bit(nt,:));
    end
    cBER_cost=sum(sum(BER));
end
return

function symbol_est=QPSK_symbol_decision(S_Est)
[Nt,K]=size(S_Est);
symbol_est=zeros(Nt,K);
for n=1:Nt
    for k=1:K
        if(real(S_Est(n,k))-imag(S_Est(n,k))>0 && real(S_Est(n,k))+imag(S_Est(n,k))>0)
            symbol_est(n,k)=1.0000+0.0000i;
        else if(real(S_Est(n,k))-imag(S_Est(n,k))>0 && real(S_Est(n,k))+imag(S_Est(n,k))<0)
                symbol_est(n,k)=0.0000-1.0000i;
            else if(real(S_Est(n,k))-imag(S_Est(n,k))<0 && real(S_Est(n,k))+imag(S_Est(n,k))>0)
                    symbol_est(n,k)=0.0000+1.0000i;
                else if(real(S_Est(n,k))-imag(S_Est(n,k))<0 && real(S_Est(n,k))+imag(S_Est(n,k))<0)
                        symbol_est(n,k)=-1.0000+0.0000i;
                    end
                end
            end
        end
    end
end       
return