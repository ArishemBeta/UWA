%Purpose: Enhanced soft-decision equalization for ovsp MIMO OFDM
%Author: Jun Tao
%
%
function [S_Est LLe_cod] = Enh_Ovsp_MIMO_OFDM_SoftEqu_qpsk(y,LLa_cod,h,K,osf,L,Nt,Nr,Nbps,SNRdB)
%function [S_Soft_Est LLe_cod] = Ovsp_MIMO_OFDM_SoftEqu_qpsk(y,LLa_cod,h,K,osf,L,Nt,Nr,Nbps,SNRdB)

%Nt: number of transmit antennas
%K: number of subcarriers

SNR= 10^(SNRdB/10);  %SNR in linear scale
if(1)
    pn= norm(h,'fro').^2/Nr/SNR/osf; %evaluation of noise power
    %pn= 1e-8;
    %pn= 1/SNR;
    SNR= 1/pn;  %(Tx Sig Pwr / Noise Pwr)
end

Q= 2^Nbps; %constellation size

QPSK_Sym_Set= [1 -i i -1]; 
QPSK_Bit_Set= [0 0;1 0;0 1;1 1];

S_Est= zeros(Nt, K);
S_Soft_Est= zeros(Nt, K);

LLa_sym= zeros(Q*Nt,K);
norm_fac= zeros(1,K);
LLa_sym_norm= zeros(Q,K);
Pa_sym= zeros(Q,K);

LLe_cod= zeros(2*Nt,K*Nbps);

for nt= 1: Nt
    %convert bit LL to symbol LL
    LLa_sym((nt-1)*Q+1,:)= LLa_cod((nt-1)*2+1,1:2:end)+LLa_cod((nt-1)*2+1,2:2:end); %1
    LLa_sym((nt-1)*Q+2,:)= LLa_cod((nt-1)*2+2,1:2:end)+LLa_cod((nt-1)*2+1,2:2:end); %-i
    LLa_sym((nt-1)*Q+3,:)= LLa_cod((nt-1)*2+1,1:2:end)+LLa_cod((nt-1)*2+2,2:2:end); %i
    LLa_sym((nt-1)*Q+4,:)= LLa_cod((nt-1)*2+2,1:2:end)+LLa_cod((nt-1)*2+2,2:2:end); %-1
    %normalize symbol a priori probability
    for m= 1: K
        norm_fac(m)= -logsum(LLa_sym((nt-1)*Q+1: nt*Q,m).');
    end
    LLa_sym_norm= LLa_sym((nt-1)*Q+1: nt*Q,:)+ ones(Q,1)*norm_fac;
    Pa_sym= exp(LLa_sym_norm);
    Pa_sym= Pa_sym./(ones(Q,1)*sum(Pa_sym,1));
    
    %calculate the mean and variance of the symbols   
    S_bar(nt,:)= QPSK_Sym_Set* Pa_sym; 
    %calculate a priori variance of Tx symbols
    S_var(nt,:)= abs(QPSK_Sym_Set).^2* Pa_sym- abs(S_bar(nt,:)).^2;
    %stack mean/variance vector of all Tx into one vector
    %ss_bar((nt-1)*Nblk+1: nt*Nblk)= s_bar(nt,:);
    %ss_var((nt-1)*Nblk+1: nt*Nblk)= s_var(nt,:);
end

%S_var(S_var<1e-8)= 1e-8; %added by Jun Tao on Oct. 18, 2010, to avoid overflow

E= 0;
for nr= 1: Nr
    Y(nr,:)= fft(y(nr,:));
    E= E+sum(abs(Y(nr,:)).^2);
end
E= E/Nr/K/osf;
sigma2= E/10^(SNRdB/10)

%obtain frequency-domain channel via FFT
if(1)
    for nr= 1: Nr
        for nt= 1: Nt
            H(nr,(nt-1)*osf*K+1: nt*osf*K)= fft(h(nr,(nt-1)*osf*L+1: nt*osf*L),osf*K);
        end
    end
else
    for nr= 1: Nr
        for nt= 1: Nt
            H(nr,(nt-1)*osf*K+1: nt*osf*K)= fft(h(nr,(nt-1)*osf*L+1: osf: nt*osf*L),osf*K);
        end
    end    
end

for k= 1: K
    Yk= [];
    Hk= [];
    
    if(0)
        temp_vec= [];
        for nr= 1: Nr
            temp_vec= [temp_vec 1 0.01*ones(1,osf-1)];
        end
        temp_matrix= diag(temp_vec);
        %c= [1, 0.0125-0.1090i,-0.0423-0.0521i,...
        %   0.0236-0.0081i,0.0597-0.0204i,-0.0111-0.0285i];         
    end
    
    %construct FD Rx vector at k
    for nosf= 1: osf
        Yk= [Yk Y(:,k+(nosf-1)*K)];
    end
    Yk= reshape(Yk.',osf*Nr,1);
    %construct Hk corresponding to Yk
    for nt= 1: Nt
        temp= [];
        for nosf= 1: osf
            temp= [temp H(:,(nt-1)*osf*K+k+(nosf-1)*K)];
        end
        Hk= [Hk reshape(temp.',osf*Nr,1)];
    end 
    size(Hk);
    
    %soft-decision equalization
    for nt= 1: Nt  
        %sigma2= 0.005;
        temp_nt= [zeros(1,nt-1) 1 zeros(1,Nt-nt)]';
        Hk_nt= Hk*temp_nt; %take the n-th column of Hk
        if(1)
            Bk= Hk*diag(S_var(:,k))*Hk'+sigma2*eye(osf*Nr);
        else
            Bk= Hk*diag(S_var(:,k))*Hk'+sigma2*temp_matrix;
        end
              
        bk_nt= Hk_nt'*inv(Bk)*Hk_nt;
        %symbol estimaion using a priori knowledge
        %S_Est(nt,k)= (Hk_nt)'*inv(Bk+(1-S_var(nt,k))*Hk_nt*Hk_nt');
        S_Est(nt,k)= (Hk_nt)'*inv(Bk)*(Yk-Hk*S_bar(:,k))/(1+(1-S_var(nt,k))*bk_nt)+...
            bk_nt*S_bar(nt,k)/(1+(1-S_var(nt,k))*bk_nt);
        %calculate conditional mean and variance of S_Est(nt,k)
        mk_nt= bk_nt/(1+(1-S_var(nt,k))*bk_nt);
        vk_nt= bk_nt*(1-S_var(nt,k)*bk_nt)/(1+(1-S_var(nt,k))*bk_nt)^2;
        vk_nt= real(vk_nt);
        
        %Calculate LLR for coded bits at time n
        rou1= real(abs(S_Est(nt,k)-1*mk_nt)^2/vk_nt); %QPSK sym= 1
        rou2= real(abs(S_Est(nt,k)-(-i)*mk_nt)^2/vk_nt); %QPSK sym= -i
        rou3= real(abs(S_Est(nt,k)-i*mk_nt)^2/vk_nt); %QPSK sym= i
        rou4= real(abs(S_Est(nt,k)-(-1)*mk_nt)^2/vk_nt); %QPSK sym= -1

        %LL for the first composing bit of QPSK symbol        
        a= logsum([-rou1 -rou3]+[LLa_sym((nt-1)*Q+1,k) LLa_sym((nt-1)*Q+3,k)]);
        b= logsum([-rou2 -rou4]+[LLa_sym((nt-1)*Q+2,k) LLa_sym((nt-1)*Q+4,k)]);    
        c= logsum([a b]); a= a- c; b= b- c;
        LLe_cod((nt-1)*2+1,(k-1)*Nbps+1)= a; %bit 1 equal to 0
        LLe_cod((nt-1)*2+2,(k-1)*Nbps+1)= b; %bit 1 equal to 1

        %LL for the second composing bit of QPSK symbol
        a= logsum([-rou1 -rou2]+[LLa_sym((nt-1)*Q+1,k) LLa_sym((nt-1)*Q+2,k)]);
        b= logsum([-rou3 -rou4]+[LLa_sym((nt-1)*Q+3,k) LLa_sym((nt-1)*Q+4,k)]);     
        c= logsum([a b]); a= a- c; b= b- c;
        LLe_cod((nt-1)*2+1,(k-1)*Nbps+2)= a; %bit 2 equal to 0
        LLe_cod((nt-1)*2+2,(k-1)*Nbps+2)= b; %bit 2 equal to 1
              
        sum1= [-rou1 -rou2 -rou3 -rou4]+ LLa_sym((nt-1)*Q+1: nt*Q,k).';
        sum1= sum1- logsum(sum1); %normalize
        Pa_sym_temp= exp(sum1.').';
        %Calculate symbol a posteriori mean and variance
        S_Soft_Est(nt,k)= QPSK_Sym_Set* Pa_sym_temp.';
        S_bar(nt,k)= S_Soft_Est(nt,k);
        S_var(nt,k)= 1- abs(S_bar(nt,k))^2;
    end    
    
%     %LS/LMMSE symbol estimation
%     if(1)
%         %S_Est(:,k)= pinv(Hk)*Yk;
%         S_Est(:,k)= inv(Hk'*Hk)*Hk'*Yk;
%     else
%         S_Est(:,k)= Hk'*inv(Hk*Hk'+pn*eye(osf*Nr))*Yk; 
%     end
end

LLe_cod= LLe_cod- LLa_cod;
    