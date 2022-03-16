function H=OMP(z,Npa,Namp,B,Ns,K,sc_idx,block_symbol,bmin,bmax,dbeta,cfo,L)
sparsity=Npa*(Namp+1);
t=zeros(1,K);
t(sc_idx)=1;
t(sc_idx+1)=1;
t(sc_idx-1)=1;
selector=zeros(K,K);
selector(logical(eye(size(selector))))=t;
zp=selector*z;
sp=selector*block_symbol;
delay_hat=[0:1/(B):L/B-1/(B)];
beta=[bmin:dbeta:bmax];
% beta=-0.0010;
Np=length(delay_hat);                             %----------Np-----------
Nb=length(beta);                                  %|
A=zeros(K,(Namp+1)*Np*Nb);                        %|
% Ablock=zeros(K/length(sc_idx),Np*Nb);           %Nb
syms f;                                           %|
G=(sin(pi*f*K/B)/(pi*f*K/B))*exp(-sqrt(-1)*pi*f*K/B);
for n=1:Namp+1
    Gn=diff(G,n-1);
    for b=1:Nb
        gamma=zeros(K,K);
        parfor m=1:K
            m
            for k=1:K
%                 (m-1)*K+k
                gamma(m,k)=subs(Gn,'f',(m-k)*B/K+(cfo-beta(b)*(13000+(m-1-K/2)*B/K))/(1+beta(b)));
            end
        end
        for p=1:Np
            (n-1)*Np*Nb+(b-1)*Np+p
            A(:,(n-1)*Np*Nb+(b-1)*Np+p)=diag(exp(-sqrt(-1)*2*pi*delay_hat(p)*B/K*[-K/2:K/2-1]))...
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
while k<=sparsity
    [Rm,ind]=max(abs(cor));
    index=[index ind];
    P=A(:,index)*inv(A(:,index)'*A(:,index))*A(:,index)';
    r=(eye(Am)-P)*zp;
    cor=A'*r;
    k=k+1;
end
xind=inv(A(:,index)'*A(:,index))*A(:,index)'*zp;
xi(index)=xind;

pos=zeros(length(index),3);
param=zeros(length(index),3);
for i=1:length(index)
    pos(i,1)=fix(index(i)/(Np*Nb))+1;                 %Namp
    pos(i,2)=fix((index(i)-Np*Nb*(pos(i,1)-1))/Np)+1;     %Nb
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

%     p1=fix(index(i)/(Np*Nb))+1;
%     p2=fix((index(i)-Np*Nb*(p1-1))/Np)+1;
%     p3=index(i)-Np*Nb*(p1-1)-Np*(p2-1);
%     if(p3==0)
%         p3=Np;
%         p2=p2-1;
%         if(p2==0);
%             p2=Nb;
%             p1=p1-1;
%         end
%     end
end

H=zeros(K,K);
for i=1:length(index)
    Gn=diff(G,param(i,1));
    parfor m=1:K
        for k=1:K
            gamma(m,k)=subs(Gn,'f',(m-k)*B/K+(cfo-param(i,2)*(13000+(m-1-K/2)*B/K))/(1+param(i,2)));
        end
    end
    H=H+xind(i)*diag(exp(-sqrt(-1)*2*pi*param(i,3)*B/K*[-K/2:K/2-1]))*gamma;
end
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
%     r = (eye(Am)-P)*b; cor=A'*r;
%     k=k+1;
% end
% %Step 5
% xind = inv(A(:,index)'*A(:,index))*A(:,index)'*b;
% x(index) = xind;
% end