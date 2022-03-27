function H=OMP(z,Npa,Namp,B,Ns,K,sc_idx,block_symbol,bmin,bmax,dbeta,cfo,L,D)
sparsity=Npa*(Namp+1);
t=zeros(1,K);
t(sc_idx)=1;
% t(sc_idx-1)=1;
% t(sc_idx+1)=1;
T=K/B;
selector=zeros(K,K);
selector(logical(eye(size(selector))))=t;
zp=selector*z;
sp=selector*block_symbol;
delay_hat=[1/B:1/(B):L/B];
beta=[bmin:dbeta:bmax];
Np=length(delay_hat);                             %----------Np-----------
Nb=length(beta);                                  %|
A=zeros(K,(Namp+1)*Np*Nb);                        %|
% Ablock=zeros(K/length(sc_idx),Np*Nb);           %Nb
syms f;                                           %|
G=(sin(pi*f*T)/(pi*f*T))*exp(-sqrt(-1)*pi*f*T);
for n=1:Namp+1
    Gn=diff(G,n-1);
    for b=1:Nb
        gamma=zeros(K,K);
        for d=1:D+1
            if(d==1)
                for m=1:K
                    gamma(m,m)=subs(Gn,'f',(cfo-beta(b)*(13000+(m-1-K/2)/T))/(1+beta(b)));
                end
            else
                for m=1:K-d+1
                    gamma(m,m+d-1)=subs(Gn,'f',(1-d)/T+(cfo-beta(b)*(13000+(m-1-K/2)/T))/(1+beta(b)));
                    gamma(m+d-1,m)=subs(Gn,'f',(d-1)/T+(cfo-beta(b)*(13000+(m-1-K/2)/T))/(1+beta(b)));
                end
            end
        end
        for p=1:Np
            (n-1)*Np*Nb+(b-1)*Np+p
            A(:,(n-1)*Np*Nb+(b-1)*Np+p)=diag(exp(-sqrt(-1)*2*pi*delay_hat(p)/T*[-K/2:K/2-1]))...
                *gamma...
                *sp;
        end
    end
end
xi=zeros((Namp+1)*Np*Nb,1);

index=[];
k=1;
[Am,An]=size(A);
r=zp;
cor=A'*r;
coro=sum(abs(cor));
while sum(abs(cor))>0.1
% while k<=sparsity*10
    if(k>Nb*Np*(Namp+1)) break; end;
    [Rm,ind]=max(abs(cor));
    index=[index ind];
%     P=A(:,index)*inv(A(:,index)'*A(:,index))*A(:,index)';
    P=A(:,index)/(A(:,index)'*A(:,index))*A(:,index)';
    r=(eye(Am)-P)*zp;
    cor=A'*r;
    k=k+1;
end
sum(abs(cor))
k
xind=(A(:,index)'*A(:,index))\A(:,index)'*zp;
xi(index)=xind;

pos=zeros(length(index),3);
param=zeros(length(index),5);
for i=1:length(index)
    pos(i,1)=fix(index(i)/(Np*Nb))+1;                       %Namp
    pos(i,2)=fix((index(i)-Np*Nb*(pos(i,1)-1))/Np)+1;       %Nb
    pos(i,3)=index(i)-Np*Nb*(pos(i,1)-1)-Np*(pos(i,2)-1);   %Np
    if(pos(i,3)==0)
        pos(i,3)=Np;
        pos(i,2)=pos(i,2)-1;
        if(pos(i,2)==0)
            pos(i,2)=Nb;
            pos(i,1)=pos(i,1)-1;
        end
    end
    param(i,1)=pos(i,1)-1;
    param(i,2)=beta(pos(i,2));
    param(i,3)=delay_hat(pos(i,3));
end
param(:,4)=index.';
param(:,5)=(1:height(param)).';
param=sortrows(param,1);
param=sortrows(param,2);

H=zeros(K,K);
Gn=diff(G,param(1,1));
gamma=zeros(K,K);
for d=1:D+1
    if(d==1)
        for m=1:K
            gamma(m,m)=subs(Gn,'f',(cfo-param(1,2)*(13000+(m-1-K/2)/T))/(1+param(1,2)));
        end
    else
        for m=1:K-d+1
            gamma(m,m+d-1)=subs(Gn,'f',(1-d)/T+(cfo-param(1,2)*(13000+(m-1-K/2)/T))/(1+param(1,2)));
            gamma(m+d-1,m)=subs(Gn,'f',(d-1)/T+(cfo-param(1,2)*(13000+(m-1-K/2)/T))/(1+param(1,2)));
        end
    end
end

for i=1:height(param)
    if((i>1)&&(param(i,1)~=param(i-1,1)))
%         Gn=diff(G,param(i,1));
        Gn=diff(Gn,1);
    end
    if((i>1)&&(param(i,2)~=param(i-1,2)))
        for d=1:D+1
            if(d==1)
                for m=1:K
                    gamma(m,m)=subs(Gn,'f',(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
                end
            else
                for m=1:K-d+1
                    gamma(m,m+d-1)=subs(Gn,'f',(1-d)/T+(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
                    gamma(m+d-1,m)=subs(Gn,'f',(d-1)/T+(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
                end
            end
        end
    end
    H=H+xind(param(i,5))*diag(exp(-sqrt(-1)*2*pi*param(i,3)/T*[-K/2:K/2-1]))*gamma;
end

% for i=1:length(index)
%     Gn=diff(G,param(i,1));
%     gamma=zeros(K,K);
%     for d=1:D+1
%         if(d==1)
%             for m=1:K
%                 gamma(m,m)=subs(Gn,'f',(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
%             end
%         else
%             for m=1:K-d+1
%                 gamma(m,m+d-1)=subs(Gn,'f',(1-d)/T+(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
%                 gamma(m+d-1,m)=subs(Gn,'f',(d-1)/T+(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
%             end
%         end
%     end
%     H=H+xind(i)*diag(exp(-sqrt(-1)*2*pi*param(i,3)/T*[-K/2:K/2-1]))*gamma;
% end
figure();
image(abs(H),'CDataMapping','scaled');
colorbar;

return

% function x = OMP(A,b,sparsity)
% %Step 1
% index = []; k = 1; [Am, An] = size(A); r = b; x=zeros(An,1);
% cor = A'*r; 
% while k <= sparsity
%     %Step 2
%     [Rm,ind] = max(abs(cor)); 
%     index = [index ind]; 
%     %Step 3
%     P = A(:,index)*inv(A(:,index)'*A(:,index))*A(:,index)';
%     r = (eye(Am)-P)*b; 
%     cor=A'*r;
%     k=k+1;
% end
% %Step 5
% xind = inv(A(:,index)'*A(:,index))*A(:,index)'*b;
% x(index) = xind;
% end


% function theta=OMP(y,A,t)
% [y_rows,y_columns] = size(y);
% [M,N]=size(A);
% theta=zeros(N,1);
% At=[];
% Pos_theta=[];
% r_n=y;
% product=A'*r_n;
% while product>0.3
%     [val,pos]=max(abs(product));
%     At=[At,A(:,pos)];
%     Pos_theta=[Pos_theta,pos];
%     A(:,pos)=zeros(M,1);
%     theta_ls=inv(At(:,1:end)'*At(:,1:end))
% end