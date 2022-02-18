%--------------------------------------------------------------------------
%  - Time Domain Channel Equalization and Symbol Detection
%  - MIMO Turbo LE Equalization
% 
%-----------------------------------------------------
%  Created by Jun Tao on Feb.15, 2009
%  Modified on March 10, 2009, avoiding overflow
%-----------------------------------------------------
%

function [x_Est LLe_cod soft_x_Equ_all]= MIMO_MMSE_LC_TurboEqu_QPSK_fn(y,LLa_cod,H,L,Nt,Nr,Nblk,Nbps,SNRdB,head,tail)
                
SNR= 10^(SNRdB/10);  %SNR in linear scale
if(1)
    pn= norm(H,'fro').^2/Nr/SNR; %evaluation of noise power
    %pn= 1e-8;
    %pn= 1/SNR;
    SNR= 1/pn;  %(Tx Sig Pwr / Noise Pwr)
end

K1=L; %K1= 2*L; %non-causal part
K2=K1; %K2= 2*L; %causal part
K= K1+K2+1;
KK= K1+K2+L;
Q= 2^Nbps;

QPSK_Sym_Set= [1 -i i -1]; 
QPSK_Bit_Set= [0 0;1 0;0 1;1 1];

LLa_sym= zeros(Q*Nt,Nblk+KK-1);
norm_fac= zeros(1,Nblk+KK-1);
LLa_sym_norm= zeros(Q,Nblk+KK-1);
Pa_sym= zeros(Q,Nblk+KK-1);


for nt= 1: Nt
    %convert bit LL to symbol LL
    LLa_sym((nt-1)*Q+1,:)= LLa_cod((nt-1)*2+1,1:2:end)+LLa_cod((nt-1)*2+1,2:2:end); %1
    LLa_sym((nt-1)*Q+2,:)= LLa_cod((nt-1)*2+2,1:2:end)+LLa_cod((nt-1)*2+1,2:2:end); %-i
    LLa_sym((nt-1)*Q+3,:)= LLa_cod((nt-1)*2+1,1:2:end)+LLa_cod((nt-1)*2+2,2:2:end); %i
    LLa_sym((nt-1)*Q+4,:)= LLa_cod((nt-1)*2+2,1:2:end)+LLa_cod((nt-1)*2+2,2:2:end); %-1
    %normalize symbol a priori probability
    for m= 1: Nblk+KK-1
        norm_fac(m)= -logsum(LLa_sym((nt-1)*Q+1: nt*Q,m).');
    end
    LLa_sym_norm= LLa_sym((nt-1)*Q+1: nt*Q,:)+ ones(Q,1)*norm_fac;
    Pa_sym= exp(LLa_sym_norm);
    Pa_sym= Pa_sym./(ones(Q,1)*sum(Pa_sym,1));
    
    %calculate the mean and variance of the symbols   
    x_bar(nt,:)= QPSK_Sym_Set* Pa_sym; 
    %calculate a priori variance of Tx symbols
    x_var(nt,:)= abs(QPSK_Sym_Set).^2* Pa_sym- abs(x_bar(nt,:)).^2;
    %stack mean/variance vector of all Tx into one vector
    %ss_bar((nt-1)*Nblk+1: nt*Nblk)= s_bar(nt,:);
    %ss_var((nt-1)*Nblk+1: nt*Nblk)= s_var(nt,:);
end

%first (K2+L-1) and last (K1) Tx symbols are assumed to be zero
if(head)
x_bar(:,1: K2+L-1)= 0;
x_var(:,1: K2+L-1)= 0;
elseif(tail)
x_bar(:,end-K1+1: end)= 0;
x_var(:,end-K1+1: end)= 0;
end
%based on mean of x, calculate mean of Rx symbol
y_bar= zeros(Nr,Nblk+KK-1);
for nr= 1: Nr
    for nt= 1: Nt
        h= H(nr,(nt-1)*L+1: nt*L);
        y_bar(nr,:)= y_bar(nr,:)+filter(h,1,x_bar(nt,:));
    end   
end
y_bar= y_bar(:,L:end); 

%based on mean of x, calculate mean of Rx symbol
y_bar= zeros(Nr,Nblk+KK-1);
for nr= 1: Nr
    for nt= 1: Nt
        h= H(nr,(nt-1)*L+1: nt*L);
        y_bar(nr,:)= y_bar(nr,:)+filter(h,1,x_bar(nt,:));
    end   
end
y_bar= y_bar(:,L:end); 
    
%---Perform MMSE-LE equalization for the whole block---
x_Est= zeros(Nt,Nblk);

%Calculate Equalizer taps for all n,k
H1= zeros(Nr,Nt*L);
for l= 1: L
    H1(:,(l-1)*Nt+1:l*Nt)= H(:, L-l+1: L: Nt*L-l+1);
end
temp= zeros(Nr,Nt);
HH= zeros(K*Nr,KK*Nt);
for k= 1: K
    HH((k-1)*Nr+1: k*Nr,:)=[repmat(temp,1,k-1) H1 repmat(temp,1,K-k)];
end

for nt= 1: Nt
    temp= [zeros(1,nt-1) 1 zeros(1,Nt-nt)];
    s(:,nt)= HH*[zeros(1,Nt*(K2+L-1)) temp zeros(1,Nt*K1)]'; %for Tx nt
end

% for nt= 1: Nt
%     v_bar(nt)= mean(x_var(nt,K2+L: Nblk+K2+L-1));
% end

V_bar= zeros(Nt*KK,Nt*KK);
if(1)
    for nt= 1: Nt
        %for n= 1: Nblk
            v_bar(nt)= mean(x_var(nt,K2+L: Nblk+K2+L-1));
        %end
        V_bar(nt:Nt:Nt*KK,nt:Nt:Nt*KK) = v_bar(nt)*eye(KK,KK);
    end
else
    for n= 1: Nblk
        V_bar= V_bar+diag(reshape(x_var(:,n: n+KK-1),KK*Nt,1));
    end
    V_bar= V_bar/Nblk;
end
V_bar;

LLe_cod= zeros(2*Nt,Nblk*Nbps);

for nt= 1: Nt
    
    ss= s(:,nt);  vv_bar= v_bar(nt); sigma2= 1/SNR;
    if(1)
        F= inv(sigma2*eye(K*Nr)+HH*V_bar*HH')*ss;
    else %reduce complexity
        HH1= HH*V_bar^(0.5);
        temp= [zeros(1,nt-1) 1 zeros(1,Nt-nt)];
        E= [zeros(1,Nt*(K2+L-1)) temp zeros(1,Nt*K1)]; %for Tx nt
        F= (E*V_bar^(-0.5)*inv(sigma2*eye(KK*Nt)+HH1'*HH1)*HH1')';
    end
    mu= F'*ss;
    mu= real(F'*ss);
    p= HH'*F;
    C= sigma2*F'*F;

    for n= 1: Nblk
        yn= reshape(y(:,n:n+K1+K2),K*Nr,1); %received vector at time n
        yn_bar= reshape(y_bar(:,n:n+K1+K2),K*Nr,1); %received mean vector at time n
        Vn= diag(reshape(x_var(:,n: n+K1+K2+L-1),KK*Nt,1));

        %perform MMSE LE equalization
        x_Est(nt,n)= F'*(yn-yn_bar)+x_bar(nt,K2+L+n-1)*mu;
        x_Est(nt,n);
        %Calculate LLR for coded bits at time n
        rou1= real(abs(x_Est(nt,n)-1*mu)^2/real(C+p'*Vn*p-x_var(nt,K2+L+n-1)*mu^2)); %QPSK sym= 1
        rou2= real(abs(x_Est(nt,n)-(-i)*mu)^2/real(C+p'*Vn*p-x_var(nt,K2+L+n-1)*mu^2)); %QPSK sym= -i
        rou3= real(abs(x_Est(nt,n)-i*mu)^2/real(C+p'*Vn*p-x_var(nt,K2+L+n-1)*mu^2)); %QPSK sym= i
        rou4= real(abs(x_Est(nt,n)-(-1)*mu)^2/real(C+p'*Vn*p-x_var(nt,K2+L+n-1)*mu^2)); %QPSK sym= -1
        
        if(0)
        C
        p'*Vn*p
        x_var(nt,K2+L+n-1)*mu^2
        C+p'*Vn*p-x_var(nt,K2+L+n-1)*mu^2
        end
        
        %LL for the first composing bit of QPSK symbol        
        a= logsum([-rou1 -rou3]+[LLa_sym((nt-1)*Q+1,K2+L-1+n) LLa_sym((nt-1)*Q+3,K2+L-1+n)]);
        b= logsum([-rou2 -rou4]+[LLa_sym((nt-1)*Q+2,K2+L-1+n) LLa_sym((nt-1)*Q+4,K2+L-1+n)]);    
        c= logsum([a b]); 
        a= a- c;
        b= b- c;
        LLe_cod((nt-1)*2+1,(n-1)*Nbps+1)= a; %bit 2 equal to 0
        LLe_cod((nt-1)*2+2,(n-1)*Nbps+1)= b; %bit 2 equal to 1
       
        %LL for the second composing bit of QPSK symbol
        a= logsum([-rou1 -rou2]+[LLa_sym((nt-1)*Q+1,K2+L-1+n) LLa_sym((nt-1)*Q+2,K2+L-1+n)]);
        b= logsum([-rou3 -rou4]+[LLa_sym((nt-1)*Q+3,K2+L-1+n) LLa_sym((nt-1)*Q+4,K2+L-1+n)]);     
        c= logsum([a b]); a= a- c; b= b- c;
        LLe_cod((nt-1)*2+1,(n-1)*Nbps+2)= a; %bit 2 equal to 0
        LLe_cod((nt-1)*2+2,(n-1)*Nbps+2)= b; %bit 2 equal to 1
    end
end

%LLe_cod= real(LLe_cod);
%subtract input a priori LLR
soft_x_Equ_all= [];
soft_x_Equ_all= SoftSymCal_QPSK(LLe_cod,Nt,Nblk);

LLe_cod= LLe_cod- LLa_cod(:,Nbps*(K2+L-1)+1: Nbps*(K2+L-1+Nblk));

return

%normalize again
for nt= 1: Nt
   norm1= zeros(1,Nblk*Nbps);
   for n= 1: Nblk*Nbps
        norm1(n)= logsum(LLe_cod((nt-1)*2+1: nt*2,n));
   end   
   LLe_cod((nt-1)*2+1: nt*2,:)= LLe_cod((nt-1)*2+1: nt*2,:)-ones(2,1)*norm1;
end

LLe_cod(LLe_cod>10)= 10;
LLe_cod(LLe_cod<-10)= -10;

return
