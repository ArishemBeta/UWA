clc;
clear all;
close all;

K=1024;
B=9765.625;
Nbps=2; %QPSK
Nt=1;
Nr=1;
Ns=4;
T214 = poly2trellis(4,[17 13]); %For convolutional encoder
tblen=12;               %For convolutional decoder
rate= 1/2;              %convolutional coding rate
sc_idx= [1:4:K];        %indices of subcarriers for channel estimation
Nbit= K*Nbps*rate;      %number of information bits in one OFDM symbol
SNRdB= 10;              %signal-to-noise ratio in dB
SNR= 10^(SNRdB/10);     %SNR in linear scale
Kg= 240;                %gap between OFDM symbol or between OFDM symbol and pilot block
L= 80;                  %length of channel
Nfrm=1;                 % frame
Nblock=2;              % block per frame
M=4;

%------------generate Tx bits--------------
P_bit=randi([0 1],1,K*Nblock*Nfrm);
% load P_bitQPSK.mat;
bit_seq=reshape(P_bit,K,length(P_bit)/K);

%------------channel coding--------------
for i=1:length(P_bit)/K
code_bitt(i,:)=convenc(bit_seq(:,i),T214);
code_bit(i,:)=randintrlv(code_bitt(i,:),0);
%code_bit=code_bitt;
end

%------------QPSK--------------
data_sym_t=qpsk(code_bit).';

%------------频域插零----------
% data_sym=data_sym_t;
data_sym=zeros(K*Ns,Nblock);
for i=1:Nblock
    for j=1:K/2
        data_sym(j,i)=data_sym_t(j,i);
    end
    for j=K/2+1:K
        data_sym(j+3*K,i)=data_sym_t(j,i);
    end
end

%------------S/P--------------
% sym_seq=reshape(data_sym,K,length(data_sym)/K);

%------------IFFT--------------
% ifft_data=ifft(sym_seq); 
ifft_data=ifft(data_sym);

%------------ZP--------------
Tx_block=[ifft_data;zeros(Kg*Ns,Nfrm*Nblock)];

%------------P/S--------------
Tx_data=reshape(Tx_block,[],1).';
% figure(1);


%------------转通带-----------
% lowpass(Tx_data,B,B*Ns);
for i=1:length(Tx_data)
    Tx_data(i)=real(Tx_data(i)*(exp(sqrt(-1)*2*pi*13000*(i-1)/(B*Ns))));%
end

%------------channel--------------
Rx_data=Tx_data;
figure();
plot(Rx_data);

%------------转基带-----------
for i=1:length(Rx_data)
    Rx_data(i)=Rx_data(i)*(exp(sqrt(-1)*2*pi*(-13000)*(i-1)/(B*Ns)));%
end
lowpass(Rx_data,B/2,B*Ns);

%------------S/P--------------
for nblk=1: 1 %Nblock
    nblk
    Noffset=(nblk-1)*(K+Kg);
    Gap=Ns*Noffset;
    block_symbol=data_sym_t(:,nblk).';
    pilot_symbol=block_symbol(sc_idx);
    block_bit=bit_seq(:,nblk).';
    pilot_bit=block_bit(sc_idx);
    Rx_block=Rx_data(:,Gap+1:Gap+round((K+Kg)*Ns));
    figure();
    plot(abs(Rx_block(1,:)));
    
    decsym=fft(Rx_block(1:Ns:K*Ns));

    
    for nt= 1: Nt
        %BER calculation for turbo detection
        Decod_InfoBit= LLR_info(nt,:)<0;
        ErrNum(nt)= sum(block_bit~=Decod_InfoBit);
        ber(nt)= ErrNum(nt)/Nbit;
        Dec_CodBit= randdeintrlv(demod_qpsk(S_Est(nt,:)),0);
        Decod_InfoBit1= vitdec(Dec_CodBit,T214,tblen,'cont','hard');
        Decod_InfoBit1= Decod_InfoBit1(tblen+1: Nbit);
        ErrNum1(nt)= sum(block_bit(1: K-tblen)~=Decod_InfoBit1);
        ErrNum2(nt)=sum(code_bitt(nblk,:)~=Dec_CodBit);
        ber1(nt)= ErrNum1(nt)/Nbit;
        ber_rec(nblk)=ber(nt);
        berraw(nt)=ErrNum2(nt)/(2*Nbit);
        ber_recraw(nblk)=berraw(nt);
    end
end