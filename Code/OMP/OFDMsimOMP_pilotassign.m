clc;
clear all;
close all;

% parpool('local',16);

MMode='QPSK';
K=1024*0.5^0;
B=9765.625*0.5^0;
ifOMP=1;
ifLS=0;

if(strcmp(MMode,'QPSK'))
    Nbps=2;
elseif(strcmp(MMode,'16QAM'))
    Nbps=4;
end
% pattern=[0,2,0,0,1,1,1,1;
%                     0,0,2,0,1,1,1,1].';
pattern=[0 2 0 0 1 1 1 0
                    0 2 0 0 1 1 1 0;
                    0 2 0 0 1 1 1 0].';
lp=height(pattern);
nKr=length(find(pattern(:,1)==0))/lp;
pKr=length(find(pattern(:,1)==2))/lp;
dKr=1-nKr-pKr;

Nt=width(pattern);
Nr=8;
Ns=16;
sc_idx=repmat(mod(find(pattern==2),lp),[1,K/lp])+repmat([0:lp:K-1],[Nt,1]);

T214 = poly2trellis(4,[17 13]);
tblen=12;               
rate= 1/2;         
Nbit_s= K*Nbps*rate;                        % bits in an OFDM symbol including null subcarriers
Nbit_d=Nbit_s*dKr;                      % useful bits in an OFDM symbol
Nbit_p=Nbit_s*pKr;

SNRdB= 10;              
SNR= 10^(SNRdB/10);

Lms=8;                                      % length of channel (ms)
L= min(fix(Lms*B/1000),width(sc_idx));     % length of channel
Kg= 3*L;                                    % guard gap
Nfrm=1;                                     % frame
Nblock=2;                                  % block per frame

%------------generate Tx bits--------------
Data_bit=randi([0 1],1,Nbit_d*Nblock*Nfrm*Nt);
Pilot_bit=randi([0 1],1,Nbit_p*Nblock*Nfrm*Nt);
% load P_bitQPSK.mat;
data_seq=reshape(Data_bit,Nbit_d,Nblock*Nfrm*Nt);
pilot_seq=reshape(Pilot_bit,Nbit_p,Nblock*Nfrm*Nt);

%------------channel coding--------------
for i=1:width(data_seq)
    code_bitt_d(i,:)=convenc(data_seq(:,i),T214);
    code_bitt_p(i,:)=convenc(pilot_seq(:,i),T214);
    code_bit_d(i,:)=randintrlv(code_bitt_d(i,:),0);
    code_bit_p(i,:)=randintrlv(code_bitt_p(i,:),0);
end

if (strcmp(MMode,'QPSK'))
    data_sym_t=qpsk(code_bit_d).';
    pilot_sym_t=qpsk(code_bit_p).';
elseif (strcmp(MMode,'16QAM'))
    for i=1:Nblock*Nt
        data_sym_t(i,:)=mqam(code_bit_d(i,:),16);
        pilot_sym_t(i,:)=mqam(code_bit_p(i,:),16);
    end
    data_sym_t=data_sym_t.';
end

%------------频域插零----------
[Tx_sym_ins, Tx_sym]=Txsymbol_arrangement(data_sym_t,pilot_sym_t,K,pattern.',Ns,Nblock);
% scatterplot(data_sym_t(:,1));

%------------IFFT--------------
ifft_data=sqrt(height(Tx_sym_ins))*ifft(Tx_sym_ins);

%------------ZP--------------
Tx_block=[ifft_data;zeros(Kg*Ns,Nfrm*Nblock*Nt)];

%------------P/S--------------
Tx_data=reshape(Tx_block,[],Nt).';

%------------转通带-----------
for i=1:width(Tx_data)
    Tx_data(:,i)=real(Tx_data(:,i)*(exp(sqrt(-1)*2*pi*13000*(i-1)/(B*Ns))));%
end
        
%------------channel--------------
Npath=8;
Rx_data=zeros(Nr,Nblock*(K+Kg)*Ns+L*Ns);
for nr=1:Nr
    for nt=1:Nt
        delay=0.5*rand(1,Npath)+[1:Npath];
        % doppler=[-0.002000,-0.002050,-0.00190,-0.001950,-0.00205000,-0.002100,-0.00200,-0.00000,...
        %     -0.00000,-0.00000,-0.00000,-0.00000,-0.00000,-0.00000,-0.00000,-0.00000];
        doppler=-3./(1500+randperm(61,Npath)-31);
        len=8;                                                      %ms
        res=0.01;                                                %ms
        UWAchannel=UWAchannel_generation(Npath,Nblock*(K+Kg)/B,1/(Ns*B),len,res,delay(1:Npath),doppler(1:Npath));
        for t=1:height(UWAchannel)
            for p=1:Npath
                if(round(t-UWAchannel(t,p*2)*Ns*B)<=0)
                    Rx_data(nr,t)=0;
                else
                    Rx_data(nr,t)=Rx_data(nr,t)+UWAchannel(t,p*2-1)*...%round(t-UWAchannel(t,p*2)*Ns)
                        Tx_data(nt,round(t-UWAchannel(t,p*2)*Ns*B));
                end
            end
        end
    end
end
Rx_datat=Rx_data;
Rx_data=awgn(Rx_data,SNRdB,'measured');
noise=Rx_data-Rx_datat;
N00=sum(abs(noise).^2)/length(noise);

%------------转基带-----------
for i=1:length(Rx_data)
    Rx_data(:,i)=Rx_data(:,i)*(exp(sqrt(-1)*2*pi*(-13000)*(i-1)/(B*Ns)));
end
Rx_data=lowpass(Rx_data.',B*0.5,B*Ns).';

%------------S/P--------------
doppler_scale=zeros(1,Nblock);
cfo=zeros(1,Nblock);
noise=[ ];

for nblk=1: 1
    nblk

    Noffset=(nblk-1)*(K+Kg)-floor(doppler_scale(max(1,nblk-1))*(K+Kg)*(nblk-1));
    Gap=Ns*Noffset;
    block_symbol=Tx_sym(:,nblk:Nblock:end);
    pilot_symbol=block_symbol(sc_idx.').';
    block_symbol=block_symbol.';
    data_bit=data_seq(:,nblk:Nblock:end).';
    Rx_block=Rx_data(:,Gap+1:Gap+round((K+L)*Ns));%/(1-0.005)                             5~0.0005    30~0.003
    
    doppler_scale(nblk)=doppler(1);
    cfo(nblk)=-0.5;%(-doppler_scale(nblk)+doppler(1))*13000+rand()-0.5;
    if (strcmp(MMode,'QPSK'))
%         [doppler_scale(nblk),cfo(nblk)]=D2SearchQPSK(doppler(1)-0.001,doppler(1)+0.001,-2,2,1,1,1,sc_idx,Nt,Nr,Ns,K,L,Rx_block,pilot_symbol,block_symbol,SNR,SNRdB,Nbps,B);
    elseif (strcmp(MMode,'16QAM'))
%         [doppler_scale(nblk),cfo(nblk)]=D2SearchQAM(-0.0031,-0.0029,0,0,1,1,1,sc_idx,Nt,Nr,Ns,K,L,Rx_block,pilot_symbol,block_symbol,SNR,SNRdB,Nbps,B);
    end
    for i=1:Nr
        Rx_block(i,:)=interp1((0:length(Rx_block(i,:))-1),Rx_block(i,:),(0:length(Rx_block(i,:))-1)/(1+doppler_scale(nblk)),'spline');
        for n=1:length(Rx_block(i,:))
            Rx_block(i,n)=Rx_block(i,n)*(exp(sqrt(-1)*2*pi*(-13000*doppler_scale(nblk))*(n-1)/(B*Ns)));
        end
    end
    for i=1:Nr
        for n=1:length(Rx_block(i,:))
            Rx_block(i,n)=Rx_block(i,n)*(exp(sqrt(-1)*2*pi*(-cfo(nblk))*(n-1)/(B*Ns)));%48828.125
        end
    end

%------------OMP------------
% N0=sum(abs(Rx_data(:,Gap+(K+L+1)*Ns:Ns:Gap+(K+Kg-10)*Ns)).^2)/((Kg-L-10));
P=sum(sum(abs(Rx_block(:,1:K*Ns)).^2))/(K*Nr*Ns);
N0=P/(SNR+1);
if(ifOMP)
    yO=Rx_block(:,1: Ns: Ns*K).';
    olaO=Rx_block(:,1+K*Ns: Ns: (K+L-1)*Ns).';
    yO(1:L-1,:)=yO(1:L-1,:)+olaO(1:L-1,:);

    z=(fft(yO)./sqrt(height(yO)));
    tic
    H=OMP_test(z,0,B,K,sc_idx.',block_symbol.',-0.00005,0.00005,0.00001,cfo(nblk),L,4,SNR);
    toc
%     S_EstO=((H'*H+N0*eye(K))\H'*z).';
    tic
    S_EstO=MIMO_LMMSE_Equalization(z,H,K,pattern,Nr,Nt,MMode,nKr,N0);
    toc
    SO=zeros(K,Nt);
    SO_data=zeros(K,Nt);
    for nt=1:Nt
        SO(:,nt)=(diag(repmat(pattern(:,nt).'>0,1,K/lp))*S_EstO(nt,:).');
        SO_data(:,nt)=(diag(repmat(pattern(:,nt).'==1,1,K/lp))*S_EstO(nt,:).');
    end
    SO(SO==0)=[];
    SO_data(SO_data==0)=[];
    SO=reshape(SO.',[],Nt).';
    SO_data=reshape(SO_data.',[],Nt).';
    scatterplot(SO(2,:));
    title('OMP',FontSize=20);

    for nt=1:Nt
        if (strcmp(MMode,'QPSK'))
            Dec_CodBitO= randdeintrlv(demod_qpsk(SO_data(nt,:)),0);
        elseif (strcmp(MMode,'16QAM'))
            Dec_CodBitO= randdeintrlv(demod_mqam(SO_data(nt,:),16),0);
        end
        Decod_InfoBitO= vitdec(Dec_CodBitO,T214,tblen,'cont','hard');
        Decod_InfoBitO= Decod_InfoBitO(tblen+1: Nbit_d);
        ErrNum1= sum(data_bit(nt,1: Nbit_d-tblen)~=Decod_InfoBitO);
        ErrNum2=sum(code_bitt_d((nt-1)*Nblock+nblk,:)~=Dec_CodBitO);
        ber_recO(nt,nblk)= ErrNum1/(Nbit_d);
        ber_recrawO(nt,nblk)=ErrNum2/(Nbit_d/rate);
    end

    BER_costO=0;
    symbol_errO=data_sym_t(:,nblk:Nblock:end).'-SO_data;
    for nt=1:Nt
        for k=1:length(symbol_errO)
            BER_costO=BER_costO+(symbol_errO(nt,k)*conj(symbol_errO(nt,k)));
        end
    end

end

%-----------LS-------------
    if(ifLS)
        y=Rx_block(:,1: Ns: Ns*K);
        ola=Rx_block(:,1+K*Ns: Ns: (K+L-1)*Ns);
        y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);
        h= FreqDomain_MIMO_ChnnEst_fn_26May11(y, block_symbol, Nt, Nr, Ns, K, L, sc_idx, SNR);
        %     figure();
        %     plot(abs(h(1,:)));

        LLa_cod= log(0.5)*ones(2*Nt,Nbps*K);
        LLR_info= zeros(Nt,Nbit_s);

        if (strcmp(MMode,'QPSK'))
            [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_qpsk_fn_28may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB);
        elseif (strcmp(MMode,'16QAM'))
            [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_16qam_fn_31may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB);
        end
        S_Est_Iter= S_Est;%((k-1)*Nt+1: k*Nt,:)

        S=(diag(repmat(pattern,1,K/lp))*S_Est.').';
        S(S==0)=[];
        S=reshape(S,Nt,[]);
        scatterplot(S(1,:));
        title('LS',FontSize=20);

        BER_cost=0;
        symbol_err=data_sym_t(:,nblk:Nblock:end).'-S;
        for nt=1:Nt
            for k=1:length(symbol_err)
                BER_cost=BER_cost+(symbol_err(nt,k)*conj(symbol_err(nt,k)));
            end
        end

        for nt=1:Nt
            if (strcmp(MMode,'QPSK'))
                Dec_CodBit= randdeintrlv(demod_qpsk(S(nt,:)),0);
            elseif (strcmp(MMode,'16QAM'))
                Dec_CodBit= randdeintrlv(demod_mqam(S(nt,:),16),0);
            end

            Decod_InfoBit1= vitdec(Dec_CodBit,T214,tblen,'cont','hard');
            Decod_InfoBit1= Decod_InfoBit1(tblen+1: Nbit_d);
            ErrNum1= sum(data_bit(nt,1: Nbit_d-tblen)~=Decod_InfoBit1);
            ErrNum2=sum(code_bitt_d((nt-1)*Nblock+nblk,:)~=Dec_CodBit);
            ber_rec(nt,nblk)=ErrNum1/(Nbit_d);
            ber_recraw(nt,nblk)=ErrNum2/(Nbit_d/rate);
        end
    end
end
% delete(gcp('nocreate'));