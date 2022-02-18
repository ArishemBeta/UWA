%Purpose: Soft-decision equalization for MIMO OFDM
%Author: Jun Tao
%Place: MST
%Date: April 29, 2011
%
%
function [S_Est S_bar_temp S_var_temp LLe_cod] = MIMO_OFDM_SoftEqu_8psk_fn_28may11(y,LLa_cod,h,K,Ns,L,Nt,Nr,Nbps,SNRdB)
%Nt: number of transmit antennas
%K: number of subcarriers

Q= 2^Nbps; %constellation size

OPSK_Sym_Set= [(-1+i)/sqrt(2), -1, i, (1+i)/sqrt(2), -i, (-1-i)/sqrt(2), (1-i)/sqrt(2), 1]; 
OPSK_Bit_Set= [0 0 0; 1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1];

S_bar= zeros(Nt, K); 
S_var= zeros(Nt, K); 
S_bar_temp= zeros(Nt, K); 
S_var_temp= zeros(Nt, K); 
LLa_sym= zeros(Q*Nt,K);
norm_fac= zeros(1,K);
LLa_sym_norm= zeros(Q,K);
Pa_sym= zeros(Q,K);

S_Est= zeros(Nt, K);
S_Soft_Est= zeros(Nt, K);
LL_cod= zeros(2*Nt,K*Nbps);

for nt= 1: Nt
    %convert bit LL to symbol LL
    LLa_sym((nt-1)*Q+1,:)= LLa_cod((nt-1)*2+1,1:3:end)+LLa_cod((nt-1)*2+1,2:3:end)+LLa_cod((nt-1)*2+1,3:3:end); %sym 1
    LLa_sym((nt-1)*Q+2,:)= LLa_cod((nt-1)*2+2,1:3:end)+LLa_cod((nt-1)*2+1,2:3:end)+LLa_cod((nt-1)*2+1,3:3:end); %sym 2
    LLa_sym((nt-1)*Q+3,:)= LLa_cod((nt-1)*2+1,1:3:end)+LLa_cod((nt-1)*2+2,2:3:end)+LLa_cod((nt-1)*2+1,3:3:end); %sym 3
    LLa_sym((nt-1)*Q+4,:)= LLa_cod((nt-1)*2+2,1:3:end)+LLa_cod((nt-1)*2+2,2:3:end)+LLa_cod((nt-1)*2+1,3:3:end); %sym 4
    LLa_sym((nt-1)*Q+5,:)= LLa_cod((nt-1)*2+1,1:3:end)+LLa_cod((nt-1)*2+1,2:3:end)+LLa_cod((nt-1)*2+2,3:3:end); %sym 5
    LLa_sym((nt-1)*Q+6,:)= LLa_cod((nt-1)*2+2,1:3:end)+LLa_cod((nt-1)*2+1,2:3:end)+LLa_cod((nt-1)*2+2,3:3:end); %sym 6
    LLa_sym((nt-1)*Q+7,:)= LLa_cod((nt-1)*2+1,1:3:end)+LLa_cod((nt-1)*2+2,2:3:end)+LLa_cod((nt-1)*2+2,3:3:end); %sym 7
    LLa_sym((nt-1)*Q+8,:)= LLa_cod((nt-1)*2+2,1:3:end)+LLa_cod((nt-1)*2+2,2:3:end)+LLa_cod((nt-1)*2+2,3:3:end); %sym 8
       
    %normalize symbol a priori probability
    for m= 1: K
        norm_fac(m)= -logsum(LLa_sym((nt-1)*Q+1: nt*Q,m).');
    end
    LLa_sym_norm= LLa_sym((nt-1)*Q+1: nt*Q,:)+ ones(Q,1)*norm_fac;
    Pa_sym= exp(LLa_sym_norm);
    Pa_sym= Pa_sym./(ones(Q,1)*sum(Pa_sym,1));
    
    %calculate the mean and variance of the symbols   
    S_bar(nt,:)= OPSK_Sym_Set* Pa_sym; 
    %calculate a priori variance of Tx symbols
    S_var(nt,:)= abs(OPSK_Sym_Set).^2* Pa_sym- abs(S_bar(nt,:)).^2;
    %stack mean/variance vector of all Tx into one vector
    %ss_bar((nt-1)*Nblk+1: nt*Nblk)= s_bar(nt,:);
    %ss_var((nt-1)*Nblk+1: nt*Nblk)= s_var(nt,:);
end

%S_var(S_var<1e-6)= 1e-6;

if(1)
E= 0;
for nr= 1: Nr
    Y(nr,:)= 1/sqrt(K)*fft(y(nr,:));
    E= E+sum(abs(y(nr,:)).^2);
end
E= E/Nr/K;
sigma2= E/10^(SNRdB/10)
end
%sigma2= 1e-5;

%obtain frequency-domain channel via FFT
if(1)
    for nr= 1: Nr
        for nt= 1: Nt
            H(nr,(nt-1)*K+1: nt*K)= fft(h(nr,(nt-1)*L+1: nt*L),K);
        end
    end
else
    for nr= 1: Nr
        for nt= 1: Nt
            H(nr,(nt-1)*K+1: nt*K)= fft(h(nr,(nt-1)*L+1: nt*L),K);
        end
    end    
end
H= H/sqrt(Ns); %taking into the oversampling factor, key!!!

A= OPSK_Sym_Set;

for k= 1: K
    Yk= [];
    Hk= [];
    temp_vec= [];
    for nr= 1: Nr
        temp_vec= [temp_vec 1];
    end
    temp_matrix= diag(temp_vec);
    
    %construct FD Rx vector at k
    Yk= Y(:,k);
    Yk= reshape(Yk.',Nr,1);
    %construct Hk corresponding to Yk
    for nt= 1: Nt
        temp= [];
        temp= [temp H(:,(nt-1)*K+k)];
        Hk= [Hk reshape(temp.',Nr,1)];
    end 
    size(Hk);
    
    %sigma2= norm(Hk,'fro').^2/Nr/10^(SNRdB/10);
    
    %soft-decision equalization
    for nt= 1: Nt  
        %sigma2= 0.005;
        temp_nt= [zeros(1,nt-1) 1 zeros(1,Nt-nt)]';
        Hk_nt= Hk*temp_nt; %take the n-th column of Hk
        if(1)
            Bk= Hk*diag(S_var(:,k))*Hk'+sigma2*eye(Nr);
        else
            Bk= Hk*diag(S_var(:,k))*Hk'+sigma2*temp_matrix;
        end
              
        bk_nt= Hk_nt'*inv(Bk)*Hk_nt;
        %symbol estimaion using a priori knowledge
        Yk-Hk*S_bar(:,k);
        %S_Est(nt,k)= (Hk_nt)'*inv(Bk+(1-S_var(nt,k))*Hk_nt*Hk_nt');
        S_Est(nt,k)= (Hk_nt)'*inv(Bk)*(Yk-Hk*S_bar(:,k))/(1+(1-S_var(nt,k))*bk_nt)+...
            bk_nt*(S_bar(nt,k))/(1+(1-S_var(nt,k))*bk_nt);
        %calculate conditional mean and variance of S_Est(nt,k)
        mk_nt= bk_nt/(1+(1-S_var(nt,k))*bk_nt);
        vk_nt= bk_nt*(1-S_var(nt,k)*bk_nt)/(1+(1-S_var(nt,k))*bk_nt)^2;
        vk_nt= real(vk_nt);
                
        %Calculate LLR for coded bits at time n
        rou1= real(abs(S_Est(nt,k)-A(1)*mk_nt)^2/vk_nt); %+1
        rou2= real(abs(S_Est(nt,k)-A(2)*mk_nt)^2/vk_nt); %-i
        rou3= real(abs(S_Est(nt,k)-A(3)*mk_nt)^2/vk_nt); %i
        rou4= real(abs(S_Est(nt,k)-A(4)*mk_nt)^2/vk_nt); %-1
        rou5= real(abs(S_Est(nt,k)-A(5)*mk_nt)^2/vk_nt); %+1
        rou6= real(abs(S_Est(nt,k)-A(6)*mk_nt)^2/vk_nt); %-i
        rou7= real(abs(S_Est(nt,k)-A(7)*mk_nt)^2/vk_nt); %i
        rou8= real(abs(S_Est(nt,k)-A(8)*mk_nt)^2/vk_nt); %-1
        
        %LL for the first composing bit of 16QAM symbol        
        a= logsum([-rou1 -rou3 -rou5 -rou7]+...
           [LLa_sym((nt-1)*Q+1,k) LLa_sym((nt-1)*Q+3,k) LLa_sym((nt-1)*Q+5,k) LLa_sym((nt-1)*Q+7,k)]);
        b= logsum([-rou2 -rou4 -rou6 -rou8]+...
           [LLa_sym((nt-1)*Q+2,k) LLa_sym((nt-1)*Q+4,k) LLa_sym((nt-1)*Q+6,k) LLa_sym((nt-1)*Q+8,k)]);     
        c= logsum([a b]); a= a- c; b= b- c;
        LL_cod((nt-1)*2+1,(k-1)*Nbps+1)= a; %bit 1 equal to 0
        LL_cod((nt-1)*2+2,(k-1)*Nbps+1)= b; %bit 1 equal to 1
        
        a= logsum([-rou1 -rou2 -rou5 -rou6]+...
           [LLa_sym((nt-1)*Q+1,k) LLa_sym((nt-1)*Q+2,k) LLa_sym((nt-1)*Q+5,k) LLa_sym((nt-1)*Q+6,k)]);
        b= logsum([-rou3 -rou4 -rou7 -rou8]+...
           [LLa_sym((nt-1)*Q+3,k) LLa_sym((nt-1)*Q+4,k) LLa_sym((nt-1)*Q+7,k) LLa_sym((nt-1)*Q+8,k)]);
        c= logsum([a b]); a= a- c; b= b- c;
        LL_cod((nt-1)*2+1,(k-1)*Nbps+2)= a; %bit 2 equal to 0
        LL_cod((nt-1)*2+2,(k-1)*Nbps+2)= b; %bit 2 equal to 1
        
        %LL for the third composing bit of 8PSK symbol
        a= logsum([-rou1 -rou2 -rou3 -rou4]+...
           [LLa_sym((nt-1)*Q+1,k) LLa_sym((nt-1)*Q+2,k) LLa_sym((nt-1)*Q+3,k) LLa_sym((nt-1)*Q+4,k)]);
        b= logsum([-rou5 -rou6 -rou7 -rou8]+...
           [LLa_sym((nt-1)*Q+5,k) LLa_sym((nt-1)*Q+6,k) LLa_sym((nt-1)*Q+7,k) LLa_sym((nt-1)*Q+8,k)]);  
        c= logsum([a b]); a= a- c; b= b- c;
        LL_cod((nt-1)*2+1,(k-1)*Nbps+3)= a; %bit 3 equal to 0
        LL_cod((nt-1)*2+2,(k-1)*Nbps+3)= b; %bit 3 equal to 1 
       
        sum1= [-rou1 -rou2 -rou3 -rou4 -rou5 -rou6 -rou7 -rou8]+ LLa_sym((nt-1)*Q+1: nt*Q,k).';
        sum1= sum1- logsum(sum1); %normalize
        Pa_sym_temp= exp(sum1.').'; 
        %Pa_sym_temp
        %size(Pa_sym_temp)
        %sum(Pa_sym_temp,1)
        %Pa_sym_temp= Pa_sym_temp./(ones(Q,1)*sum(Pa_sym_temp,1));
        %Calculate a posteriori mean
        S_Soft_Est(nt,k)= OPSK_Sym_Set* Pa_sym_temp.';       
        %calculate a posterioir knowledge, key!!!
        S_bar_temp(nt,k)= S_Soft_Est(nt,k);
        S_var_temp(nt,k)= 1- abs(S_bar_temp(nt,k))^2;
    end    
    
    %LS/LMMSE symbol estimation
%      if(1)
%          %S_Est(:,k)= pinv(Hk)*Yk;
%          S_Est(:,k)= inv(Hk'*Hk)*Hk'*Yk;
%      else
%          S_Est(:,k)= Hk'*inv(Hk*Hk'+pn*eye(Nr))*Yk; 
%      end
end

LLe_cod= real(LL_cod- LLa_cod);
    