function H=OMP_test(z,Namp,B,K,sc_idx,block_symbol,bmin,bmax,dbeta,cfo,L,D,SNR)
% sparsity=Npa*(Namp+1);
t=ones(1,K);
% t=zeros(1,K);
% t(sc_idx)=1;
% t(max(sc_idx-1,1))=1;
% t(sc_idx+1)=1;

Nt=width(block_symbol);
Nr=width(z);
H=zeros(Nr*K,Nt*K);
gamma=zeros(K,K);

T=K/B;
deltaf=1/T;
selector=zeros(K,K);
selector(logical(eye(size(selector))))=t;
zp=selector*z;
sp=selector*block_symbol;
delay_hat=1/B:1/(B):L/B;
beta=bmin:dbeta:bmax;
Np=length(delay_hat);                             %----------Np-----------
Nb=length(beta);                                        %|
A=zeros(K,Nt*(Namp+1)*Np*Nb);       %|Nb
syms f;                                                             %|
G=(sin(pi*f*T)/(pi*f*T))*exp(-sqrt(-1)*pi*f*T);
gamma_rec=zeros((Namp+1)*Nb,(2*D+1)*K);

%%%               gf    sparse    rec    par
for k=1:(Namp+1)*Nb
    n=floor((k-1)/Nb);
    b=mod(k-1,Nb)+1;
    Gn=diff(G,n);
    gf=matlabFunction(Gn);
    for s=1:(1+2*D)*K
        t=floor((s-1)/K);
        d=sign(mod(t,2)-0.5)*round(t/2);
        m=mod(s-1,K)+1;
        gamma_rec(k,s)=gf(-d*deltaf+(cfo-beta(b)*(13000+(m+abs(d)*logical(sign(d)-1)-1-K/2)*deltaf))/(1+beta(b)));
    end
end

for nr=1:Nr
    gamma_pos=zeros((2*D+1)*K-D*(D+1),2);
    gamma_val=zeros((2*D+1)*K-D*(D+1),1);
    for n=0:Namp
        for b=1:Nb
            for d=0:D
                if(d==0)
                    gamma_pos(1:K,:)=[1:K;1:K].';
                    gamma_val(1:K)=gamma_rec(Nb*n+b,1:K);
                else
                    gamma_pos((2*d-1)*K-(d-1)*d+1:2*d*K-d*d,:)=[1:K-d;d+1:K].';
                    gamma_pos(2*d*K-d*d+1:(2*d+1)*K-(d+1)*d,:)=[d+1:K;1:K-d].';
                    gamma_val((2*d-1)*K-(d-1)*d+1:2*d*K-d*d)=gamma_rec(Nb*n+b,(2*d-1)*K+1:2*d*K-d);
                    gamma_val(2*d*K-d*d+1:(2*d+1)*K-(d+1)*d)=gamma_rec(Nb*n+b,2*d*K+1:(2*d+1)*K-d);
                end
            end
            gamma=sparse(gamma_pos(:,1),gamma_pos(:,2),gamma_val);
            for nt=1:Nt
                gamma_sp=gamma*sp(:,nt);
                for p=1:Np
                    A(:,(nt-1)*(Namp+1)*Np*Nb+n*Np*Nb+(b-1)*Np+p)=sparse(1:K,1:K,exp(-sqrt(-1)*2*pi*delay_hat(p)*deltaf*[-K/2:K/2-1]))...
                        *gamma_sp;
                end
            end
        end
    end


    xi=zeros(Nt*(Namp+1)*Np*Nb,1);
    index=[];
    k=1;
    [Am,~]=size(A);
    E=eye(Am);
    Z=zp(:,nr);
    r=Z;
    cor=A'*r;
    cor_last=cor;

    % while norm
    %     if(k>Nb*Np*(Namp+1)) break; end;
    %     [Rm,ind]=max(abs(cor));
    %     index=[index ind];
    %     P=A(:,index)/(A(:,index)'*A(:,index))*A(:,index)';
    %     r=(eye(Am)-P)*zp;
    %     cor=A'*r;
    %     k=k+1;
    % end

    threshold=1/(1+SNR)*norm(cor,'fro');
    while norm(cor,'fro')>threshold
        if(k>Nb*Np*(Namp+1)) break; end;
        [~,ind]=max(abs(cor));
        index=[index ind];
        Atemp=A(:,index);
        P=Atemp/(Atemp'*Atemp)*Atemp';
%         P=A(:,index)/(A(:,index)'*A(:,index))*A(:,index)';
        r=(E-P)*Z;
        cor=A'*r;
        k=k+1;
    end

    xind=(A(:,index)'*A(:,index))\A(:,index)'*Z;
    xi(index)=xind;

    pos=zeros(length(index),4);
    param=zeros(length(index),6);
    for i=1:length(index)
        pos(i,1)=fix(index(i)/((Namp+1)*Np*Nb))+1;      %Nt
        pos(i,2)=fix((index(i)-(Namp+1)*Np*Nb*(pos(i,1)-1))/(Np*Nb))+1;                       %Namp
        pos(i,3)=fix((index(i)-(Namp+1)*Np*Nb*(pos(i,1)-1)-Np*Nb*(pos(i,2)-1))/Np)+1;       %Nb
        pos(i,4)=index(i)-(Namp+1)*Np*Nb*(pos(i,1)-1)-Np*Nb*(pos(i,2)-1)-Np*(pos(i,3)-1);   %Np
        if(pos(i,4)==0)
            pos(i,4)=Np;
            pos(i,3)=pos(i,3)-1;
            if(pos(i,3)==0)
                pos(i,3)=Nb;
                pos(i,2)=pos(i,2)-1;
                if(pos(i,2)==0)
                    pos(i,2)=Namp+1;
                    pos(i,1)=pos(i,1)-1;
                end
            end
        end
        param(i,1)=param(i,1);
        param(i,2)=pos(i,2)-1;
        param(i,3)=beta(pos(i,3));
        param(i,4)=delay_hat(pos(i,4));
    end
    param(:,5)=index.';
    param(:,6)=(1:height(param)).';
    param=sortrows(param,1);
    param=sortrows(param,2);
    param=sortrows(param,3);

    Htemp=zeros(K,K);
    i=0;
%     while i<=height(param)
%         i=i+1;
%         if(i==1 || param(i,2)~=param(i-1,2) || param(i,3)~=param(i-1,3))
%             for d=0:D
%                 if(d==0)
%                     gamma_val(1:K)=gamma_rec(Nb*param(i,2)+round((param(i,3)-bmin)/dbeta+1),1:K);
%                 else
%                     gamma_val((2*d-1)*K-(d-1)*d+1:2*d*K-d*d)=gamma_rec(Nb*param(i,2)+round((param(i,3)-bmin)/dbeta+1),(2*d-1)*K+1:2*d*K-d);
%                     gamma_val(2*d*K-d*d+1:(2*d+1)*K-(d+1)*d)=gamma_rec(Nb*param(i,2)+round((param(i,3)-bmin)/dbeta+1),2*d*K+1:(2*d+1)*K-d);
%                 end
%             end
%             gamma=sparse(gamma_pos(:,1),gamma_pos(:,2),gamma_val);
%         end
%         Htemp=Htemp+xind(param(i,6))...
%                             *sparse(1:K,1:K,exp(-sqrt(-1)*2*pi*param(i,4)/T*[-K/2:K/2-1]))...
%                             *gamma;
%     end

    for i=1:height(param)
        if(i==1 || param(i,2)~=param(i-1,2) || param(i,3)~=param(i-1,3))
            for d=0:D
                if(d==0)
                    gamma_val(1:K)=gamma_rec(Nb*param(i,2)+round((param(i,3)-bmin)/dbeta+1),1:K);
                else
                    gamma_val((2*d-1)*K-(d-1)*d+1:2*d*K-d*d)=gamma_rec(Nb*param(i,2)+round((param(i,3)-bmin)/dbeta+1),(2*d-1)*K+1:2*d*K-d);
                    gamma_val(2*d*K-d*d+1:(2*d+1)*K-(d+1)*d)=gamma_rec(Nb*param(i,2)+round((param(i,3)-bmin)/dbeta+1),2*d*K+1:(2*d+1)*K-d);
                end
            end
            gamma=sparse(gamma_pos(:,1),gamma_pos(:,2),gamma_val);
        end
        Htemp=Htemp+xind(param(i,6))...
                            *sparse(1:K,1:K,exp(-sqrt(-1)*2*pi*param(i,4)/T*[-K/2:K/2-1]))...
                            *gamma;
        if(i==height(param) || param(i,1)~=param(i+1,1))
            H((nr-1)*K+1:nr*K,param(i,1)*K+1:(param(i,1)+1)*K)=Htemp;
            Htemp=zeros(K,K);
        end
    end
end

% figure();
% image(abs(H),'CDataMapping','scaled');
% colorbar;

return
