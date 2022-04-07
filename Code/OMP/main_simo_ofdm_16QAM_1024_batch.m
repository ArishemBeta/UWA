 %Purpose: Process SPACE08 OFDM data: 2Tx,K=1024,QPSK
%Author: Jun Tao


clc
clear all
close all
addpath DopplerScaleEstimation/Code/RefSig
addpath DopplerScaleEstimation/Code/SPACE08_Moving
addpath DopplerScaleEstimation/Code/Turbo_OFDM_QPSK_2Tx_28May11

% addpath ../RefSig
% addpath ../SPACE08_Moving
% addpath ../Turbo_OFDM_QPSK_2Tx_28May11
% addpath ..\QAM16

Fs=43669.5;  %39.0625 kHz, Sampling Frequency
Fc=13e3;     %carrier frequency  (Hz)
Fb=9765.625;     %9.7656 kHz, Signal Bandwidth
Ns=5;
Nbps= 4;          %16QAM
Nt= 1;            %number of transducers
T214 = poly2trellis(4,[17 13]); %For convolutional encoder
tblen=12;               %For convolutional decoder
rate= 1/2;              %convolutional coding rate
K= 1024;                %number of subcarriers
Npilot= 90;            %number of pilot symbols for each OFDM symbol
Kg= 120;                %gap between OFDM symbol or between OFDM symbol and pilot block
sc_idx= [1:4:K];        %indices of subcarriers for channel estimation
Nbit= K*Nbps*rate;      %number of information bits in one OFDM symbol
Nrepeat= 12;             %number of OFDM symbols in each packet
SNRdB= 30;              %signal-to-noise ratio in dB
SNR= 10^(SNRdB/10);     %SNR in linear scale
L= 80;                 %length of channel

% Chnn_idx=[1];                              %01 phone
% Chnn_idx=[1 2];                           %02 phones
% Chnn_idx=[1 3 5 7];                      %04 phones
% Chnn_idx=[1 2 3 6 7 8];                  %06 phones
Chnn_idx=[1 2 3 4 5 6 7 8];             %08 phones
Nr= length(Chnn_idx);   %# of receive hydrophones

LL_min= -1e5;
pilot_LLR= 0; %1: using pilot LLR; 0: not using pilot LLR
Enh= 0; %1: enhanced equalization; 0: original equalization
Niter= 1;

%------------load Tx symbols/bits--------------
load One_Tx_Ant1_OFDM_16QAM_K1024.mat
QAM16_TxSym= TxSym_1Tx;
InfoBit=TxInfoBit_1Tx;
CodBit=TxCodIntrlvBit_1Tx;
clear One_Tx_Ant1_OFDM_16QAM_K1024;
load pilot_1Tx_OFDM.mat;

doppler_scale=zeros(1,12);
cfo=zeros(1,12);
err_block=zeros(1,12);
shift=zeros(1,12);

for packet_idx= [1] %[1 2 3 4 5 6 8 11 14 20 23 26 29 32]
    
    BB_name= ['RX_1Tx_OFDM_1024_BB_data_RemusRun01S05F10']
    eval(['load ' BB_name]);
    for i=1:Nr
        RX_data_n(i,:)=interp1([0:length(RX_data(Chnn_idx(i),:))-1],RX_data(Chnn_idx(i),:),[0:855/956:length(RX_data(Chnn_idx(i),:))-1],'spline');
        %  plot(abs(RX_data_nn(1,:)));
        RX_data_nn(i,:)=interp1((0:length(RX_data_n(i,:))-1),RX_data_n(i,:),(0:length(RX_data_n(i,:))-1)/(1+(Doppler(1)-0.5)/13000),'spline');
        for n=1:length(RX_data_nn(i,:))
            RX_data_nn(i,n)=RX_data_nn(i,n)*(exp(sqrt(-1)*2*pi*(-13000*(Doppler(1)-0.5)/13000)*(n-1)/48828.125));
        end
% RX_data_nn=RX_data_n;
    end
    N_shift_symbol= 6;
    N_shift_sample= 0;
    
    ErrNum= zeros(Nt,Niter);
    ErrNum1= zeros(Nt,Niter);
    BER= zeros(Nt,Niter);
    BER1= zeros(Nt,Niter);
    
    ErrNum_All= zeros(Niter,Nrepeat);
    ErrNum1_All= zeros(Niter,Nrepeat);
    BER_All= zeros(Niter,Nrepeat);
    BER1_All= zeros(Niter,Nrepeat);
    
    ErrNum_All_Tx= zeros(Niter*Nt,Nrepeat);
    
    S_Est_Iter= zeros(Nt*Niter,K);
    
    for nblk= [1:1]%: 12]%Nrepeat]
        
        nblk
        Noffset= 1300*Nt+511+189+12*2*1354;
        Gap=Ns*(Noffset+N_shift_symbol)+N_shift_sample;
        Gap=fix(Gap);
        block_symbol=QAM16_TxSym(nblk,:);
        pilot_symbol=block_symbol(sc_idx);
        block_bit=InfoBit(nblk,:);
        pilot_bit=block_bit(sc_idx);
        %-------initial channel estimation with pilots-------

        shift=fliplr([-11:0]);
        Gap= Gap+(Npilot+Kg)*Ns;
        if(nblk>1)
            Gap=Gap+shift(nblk);
            for i=2:nblk
                Gap= Gap+fix(Ns*1354/(1+doppler_scale(i-1)));
            end
            RX_block=RX_data_nn(:,Gap:fix(Gap+(K*Ns+(L)*Ns)/(1+doppler_scale(nblk-1)))-1);
        else
            RX_block=RX_data_nn(:,Gap:Gap+K*Ns+(L)*Ns-1);
        end

        for i=1:Nr
            RX_block(i,:)=interp1((0:length(RX_block(i,:))-1),RX_block(i,:),(0:length(RX_block(i,:))-1)/(1+doppler_scale(nblk)),'spline');
            for n=1:length(RX_block(i,:))
                RX_block(i,n)=RX_block(i,n)*(exp(sqrt(-1)*2*pi*(-13000*doppler_scale(nblk))*(n-1)/48828.125));
            end
            for n=1:length(RX_block(i,:))
                RX_block(i,n)=RX_block(i,n)*(exp(sqrt(-1)*2*pi*(-cfo(nblk))*(n-1)/48828.125));
            end
        end


%-------------OMP--------------
        ifOMP=1;
        if(ifOMP)
            yO=RX_block(:,1: Ns: Ns*K).';
            olaO=RX_block(:,1+K*Ns: Ns: (K+L-1)*Ns).';
            yO(1:L-1,:)=(yO(1:L-1,:)+olaO(1:L-1,:));

            z=(fft(yO)./sqrt(height(yO)));
            H=zeros(K*Nr,K);
            tic
%             for i=1:Nr
%                 H((i-1)*K+1:i*K,:)=OMP_MIMO(z(:,i),1,Fb,K,sc_idx,block_symbol.',-0.00010,0.00010,0.00001,-1,L,9,SNR);
%             end
            H=OMP_test(z,1,Fb,K,sc_idx,block_symbol.',-0.00010,0.00010,0.00001,-1,L,9,SNR);
            toc
%-------------均衡--------------
            %  S_EstO=((H'*H+N0*eye(K))\H'*z).';
            P=sum(sum(abs(RX_block).^2))/((K+L)*Nr*Ns);
            N0=P/(1+SNR);
            tic
            S_EstO=MIMO_LMMSE_Equalization(reshape(z,1,K*Nr),H,K,[1,1,1,1],Nr,Nt,'16QAM',0,N0);
            toc
            SO=S_EstO.';

            Dec_CodBitO= randdeintrlv(demod_mqam(SO,16),0);

            Decod_InfoBitO= vitdec(Dec_CodBitO,T214,tblen,'cont','hard');
            Decod_InfoBitO= Decod_InfoBitO(tblen+1: Nbit);
            ErrNum1= sum(block_bit(1: Nbit-tblen)~=Decod_InfoBitO);
            ErrNum2=sum(CodBit(nblk,:)~=demod_mqam(SO,16));
            ber_recO(nblk)= ErrNum1/(Nbit);
            ber_recrawO(nblk)=ErrNum2/(Nbit/rate);
            scatterplot(SO);
            title('OMP',FontSize=20);

            BER_costO=0;
            symbol_errO=block_symbol-SO;
            for nt=1:Nt
                for k=1:length(symbol_errO)
                    BER_costO=BER_costO+(symbol_errO(nt,k)*conj(symbol_errO(nt,k)));
                end
            end
        end



%----------LS----------
        y=RX_block(:,1: Ns: K*Ns);
        ola=RX_block(:,1+K*Ns: Ns: K*Ns+(L-1)*Ns);
        y(:,1:L-1)=y(:,1:L-1)+ola(:,1:L-1);

        FDSym= QAM16_TxSym(nblk,:);
        h= FreqDomain_MIMO_ChnnEst_fn_26May11(y, FDSym, Nt, Nr, Ns, K, L, sc_idx, SNR);
        %check CIR estimation
        if(0)
            PlotChannel(h,Nt,Nr,L);
            return
        end
        figure();
        plot(abs(h(1,:)));
        %starting MIMO equalization for OFDM symbol
        Yblk=y;

        %initial LL for coded bits and information bits
        LLa_cod= log(0.5)*ones(2*Nt,Nbps*K);
        LLR_info= zeros(Nt,Nbit);
        if(pilot_LLR)
            %replace estimated LLs of pilot bits with their true LLs.
            LLa_p= zeros(2*Nt,Nbps*length(sc_idx));
            for nt= 1: Nt
                pilot_bit= demod_mqam(QAM16_TxSym(nt,(nblk-1)*K+sc_idx),16);
                pilot_bit(pilot_bit==0)=0;
                pilot_bit(pilot_bit==1)=LL_min;
                LLa_p((nt-1)*2+1,:)= pilot_bit; %LL for bit equal to 0
                LLa_p((nt-1)*2+2,:)= LL_min-pilot_bit; %LL for bit equal to 1
            end
            %normalization
            for nt= 1: Nt
                for n= 1: length(sc_idx)*Nbps
                    norm1(n)= logsum(LLa_p((nt-1)*2+1: nt*2,n));
                end
                LLa_p((nt-1)*2+1: nt*2,:)= LLa_p((nt-1)*2+1: nt*2,:)-ones(2,1)*norm1;
            end
        end
        for k= 1: Niter
            %             k
            if(pilot_LLR)
                %replace estimated LLs of pilot bits with their true LLs.
                for nbps= 1: Nbps
                    LLa_cod(:,(sc_idx-1)*Nbps+nbps)= LLa_p(:,nbps:Nbps:end);
                end
            end
            if(~Enh) %original soft equalization
                [S_Est LLe_cod]= MIMO_OFDM_SoftEqu_16qam_fn_31may11(Yblk,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB);
            else  %enhanced soft equalization
                [S_Est LLe_cod]= MIMO_OFDM_EnhSoftEqu_16qam_fn_04June11(Yblk,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB);
            end
            S_Est_Iter((k-1)*Nt+1: k*Nt,:)= S_Est;

            scatterplot(S_Est);
            title('LS',FontSize=20);

            if(pilot_LLR)
                %replace estimated LLs of pilot bits with their true LLs.
                for nbps= 1: Nbps
                    LLe_cod(:,(sc_idx-1)*Nbps+nbps)= LLa_p(:,nbps:Nbps:end);
                end
            end
            for nt= 1: Nt
                LLe_cod((nt-1)*2+1,:)= randdeintrlv(LLe_cod((nt-1)*2+1,:), 0); %interleaving before delivering to equalizer
                LLe_cod((nt-1)*2+2,:)= randdeintrlv(LLe_cod((nt-1)*2+2,:), 0);

                %perform MAP convolutional decoding with soft information fed back
                [LLR_info(nt,:) LLe_cod((nt-1)*2+1: nt*2,:)] = MAPConvDecoder_R1(LLe_cod((nt-1)*2+1: nt*2,:),Nbit);

                %LLRe_cod_temp1(nt,:)= LLRe_cod(nt,:);
                LLe_cod((nt-1)*2+1,:)= randintrlv(LLe_cod((nt-1)*2+1,:), 0); %interleaving before delivering to equalizer
                LLe_cod((nt-1)*2+2,:)= randintrlv(LLe_cod((nt-1)*2+2,:), 0);
                LLa_cod((nt-1)*2+1: nt*2,:)= LLe_cod((nt-1)*2+1: nt*2,:);
            end

            if(1)
                for nt= 1: Nt
                    %BER calculation for turbo detection
                    Decod_InfoBit= LLR_info(nt,:)<0;
                    ErrNum(nt,k)= sum(block_bit~=Decod_InfoBit);
                    BER(nt,k)= ErrNum(nt,k)/Nbit;

                    Dec_CodBit= randdeintrlv(demod_mqam(S_Est(nt,:),16),0);
                    ErrNumraw(nt,k)=sum(TxCodIntrlvBit_1Tx(nblk,1: Nbps*K)~=demod_mqam(S_Est(nt,:),16));
                    Decod_InfoBit1= vitdec(Dec_CodBit,T214,tblen,'cont','hard');
                    Decod_InfoBit1= Decod_InfoBit1(tblen+1: Nbit);
                    ErrNum1(nt,k)= sum(block_bit(1: Nbit-tblen)~=Decod_InfoBit1);
                    BER1(nt,k)= ErrNum1(nt,k)/Nbit;
                    ber_rec(nblk)=ErrNum1(nt,k)/Nbit;
                    ber_recraw(nblk)=ErrNumraw(nt,k)/(2*Nbit);
                end
            end
        end %for loop of k
        ErrNum_All(:,nblk)= sum(ErrNum);
        ErrNum1_All(:,nblk)= sum(ErrNum1);
        %         ErrNumraw_All(:,nblk)= sum(ErrNumraw);
        %         for nt= 1: Nt
        %             ErrNum_All_Tx((nt-1)*Niter+1: nt*Niter,nblk)= ErrNum1(nt,:);
        %         end
    end
    
    BER_All= ErrNum1_All/Nbit
    ErrNum_All_Tx=[ErrNum_All_Tx sum(ErrNum_All_Tx,2)]
    %name= ['ErrNum_QPSK_K1024_' packet_name(packet_idx,:)];
    %eval(['save Results\' name '.mat ErrNum_All_Tx']);
    %aa= ErrNum_All_Tx(1:Niter,end)+ErrNum_All_Tx(Niter+1: Nt*Niter,end)
    %packet_error= [packet_error aa];
end

return
