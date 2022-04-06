function H=OMP_MIMO(z,Namp,B,K,sc_idx,block_symbol,bmin,bmax,dbeta,cfo,L,D,SNR)
% sparsity=Npa*(Namp+1);
t=ones(1,K);
% t=zeros(1,K);
% t(sc_idx)=1;
% t(max(sc_idx-1,1))=1;
% t(sc_idx+1)=1;

Nt=width(block_symbol);
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

gamma_pos=zeros((2*D+1)*K-D*(D+1),2);
gamma_val=zeros((2*D+1)*K-D*(D+1),1);
for nt=1:Nt
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
            for p=1:Np
                A(:,(nt-1)*(Namp+1)*Np*Nb+n*Np*Nb+(b-1)*Np+p)=sparse(1:K,1:K,exp(-sqrt(-1)*2*pi*delay_hat(p)*deltaf*[-K/2:K/2-1]))...
                    *gamma...
                    *sp(:,nt);
            end
        end
    end
end

% for n=0:Namp
%     for b=1:Nb
%         for d=0:D
%             if(d==0)
%                 for m=1:K
%                     gamma(m,m)=gamma_rec(Nb*n+b,m);
%                 end
%             else
%                 for m=1:K-d
%                     gamma(m,m+d)=gamma_rec(Nb*n+b,(2*d-1)*K+m);
%                     gamma(m+d,m)=gamma_rec(Nb*n+b,2*d*K+m);
% %                     gamma(m,m+d)=gamma_rec(Nb*n+b,(2*d-1)*K-d*(d-1)+m);
% %                     gamma(m+d,m)=gamma_rec(Nb*n+b,2*d*K-d^2+m);
%                 end
%             end
%         end
%         for p=1:Np
%             A(:,n*Np*Nb+(b-1)*Np+p)=sparse(1:K,1:K,exp(-sqrt(-1)*2*pi*delay_hat(p)*deltaf*[-K/2:K/2-1]))...
%                 *sparse(gamma)...
%                 *sp;
%         end
%     end
% end
% toc

%%%               rec

% for n=0:Namp
%     Gn=diff(G,n);
%     for b=1:Nb
%         b
%         for d=0:D
%             if(d==0)
%                 for m=1:K
%                     gamma_rec(Nb*n+b,m)=subs(Gn,'f',(cfo-beta(b)*(13000+(m-1-K/2)/T))/(1+beta(b)));
%                     gamma(m,m)= gamma_rec(Nb*n+b,m);
%                 end
%             else
%                 for m=1:K-d
%                     gamma_rec(Nb*n+b,(2*d-1)*K-d*(d-1)+m)=subs(Gn,'f',(-d)/T+(cfo-beta(b)*(13000+(m-1-K/2)/T))/(1+beta(b)));
%                     gamma_rec(Nb*n+b,2*d*K-d^2+m)=subs(Gn,'f',d/T+(cfo-beta(b)*(13000+(m+d-1-K/2)/T))/(1+beta(b)));
%                     gamma(m,m+d)=gamma_rec(Nb*n+b,(2*d-1)*K-d*(d-1)+m);
%                     gamma(m+d,m)=gamma_rec(Nb*n+b,2*d*K-d^2+m);
%                 end
%             end
%         end
%         for p=1:Np
% %             n*Np*Nb+(b-1)*Np+p
%             A(:,n*Np*Nb+(b-1)*Np+p)=diag(exp(-sqrt(-1)*2*pi*delay_hat(p)/T*[-K/2:K/2-1]))...
%                 *gamma...
%                 *sp;
%         end
%     end
% end

% % %              gf

% % % for n=1:Namp+1
% % %     Gn=diff(G,n-1);
% % %     gf=matlabFunction(Gn);
% % %     for b=1:Nb
% % %         gamma=zeros(K,K);
% % %         for d=1:D+1
% % %             if(d==1)
% % %                 for m=1:K
% % %                     gamma(m,m)=gf((cfo-beta(b)*(13000+(m-1-K/2)/T))/(1+beta(b)));
% % %                 end
% % %             else
% % %                 for m=1:K-d+1
% % %                     gamma(m,m+d-1)=gf((1-d)/T+(cfo-beta(b)*(13000+(m-1-K/2)/T))/(1+beta(b)));
% % %                     gamma(m+d-1,m)=gf((d-1)/T+(cfo-beta(b)*(13000+(m-1-K/2)/T))/(1+beta(b)));
% % %                 end
% % %             end
% % %         end
% % %         for p=1:Np
% % %             (n-1)*Np*Nb+(b-1)*Np+p
% % %             A(:,(n-1)*Np*Nb+(b-1)*Np+p)=diag(exp(-sqrt(-1)*2*pi*delay_hat(p)/T*[-K/2:K/2-1]))...
% % %                 *gamma...  
% % %                 *sp;
% % %         end
% % %     end
% % % end

xi=zeros(Nt*(Namp+1)*Np*Nb,1);
index=[];
k=1;
[Am,An]=size(A);
r=zp;
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
    [Rm,ind]=max(abs(cor));
    index=[index ind];
    P=A(:,index)/(A(:,index)'*A(:,index))*A(:,index)';
    r=(eye(Am)-P)*zp;
    cor=A'*r;
    k=k+1;
end

xind=(A(:,index)'*A(:,index))\A(:,index)'*zp;
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

H=zeros(Nt*K,K);
gamma=zeros(K,K);

for i=1:height(param)
    for d=0:D
        if(d==0)
            for m=1:K
                gamma(m,m)=gamma_rec(Nb*param(i,2)+round((param(i,3)-bmin)/dbeta+1),m);
            end
        else
            for m=1:K-d
                gamma(m,m+d)=gamma_rec(Nb*param(i,2)+round((param(i,3)-bmin)/dbeta+1),(2*d-1)*K+m);
                gamma(m+d,m)=gamma_rec(Nb*param(i,2)+round((param(i,3)-bmin)/dbeta+1),2*d*K+m);
            end
        end
    end
    H(param(i,1)*K+1:(param(i,1)+1)*K,:)=H(param(i,1)*K+1:(param(i,1)+1)*K,:)...
                                                                                    +xind(param(i,6))...
                                                                                    *sparse(1:K,1:K,exp(-sqrt(-1)*2*pi*param(i,4)/T*[-K/2:K/2-1]))...
                                                                                    *sparse(gamma);
end

% % % Gn=diff(G,param(1,1));
% % % 
% % % for d=1:D+1
% % %     if(d==1)
% % %         for m=1:K
% % %             gamma(m,m)=subs(Gn,'f',(cfo-param(1,2)*(13000+(m-1-K/2)/T))/(1+param(1,2)));
% % %         end
% % %     else
% % %         for m=1:K-d+1
% % %             gamma(m,m+d-1)=subs(Gn,'f',(1-d)/T+(cfo-param(1,2)*(13000+(m-1-K/2)/T))/(1+param(1,2)));
% % %             gamma(m+d-1,m)=subs(Gn,'f',(d-1)/T+(cfo-param(1,2)*(13000+(m-1-K/2)/T))/(1+param(1,2)));
% % %         end
% % %     end
% % % end
% % % 
% % % for i=1:height(param)
% % %     if((i>1)&&(param(i,1)~=param(i-1,1)))
% % % %         Gn=diff(G,param(i,1));
% % %         Gn=diff(Gn,1);
% % %     end
% % %     if((i>1)&&(param(i,2)~=param(i-1,2)))
% % %         for d=1:D+1
% % %             if(d==1)
% % %                 for m=1:K
% % %                     gamma(m,m)=subs(Gn,'f',(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
% % %                 end
% % %             else
% % %                 for m=1:K-d+1
% % %                     gamma(m,m+d-1)=subs(Gn,'f',(1-d)/T+(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
% % %                     gamma(m+d-1,m)=subs(Gn,'f',(d-1)/T+(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
% % %                 end
% % %             end
% % %         end
% % %     end
% % %     H=H+xind(param(i,5))*diag(exp(-sqrt(-1)*2*pi*param(i,3)/T*[-K/2:K/2-1]))*gamma;
% % % end

% % % % for i=1:length(index)
% % % %     Gn=diff(G,param(i,1));
% % % %     gf=matlabFunction(Gn);
% % % %     gamma=zeros(K,K);
% % % %     for d=1:D+1
% % % %         if(d==1)
% % % %             for m=1:K
% % % %                 gamma(m,m)=gf((cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
% % % % %                 gamma(m,m)=subs(Gn,'f',(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
% % % %             end
% % % %         else
% % % %             for m=1:K-d+1
% % % %                 gamma(m,m+d-1)=gf((1-d)/T+(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
% % % %                 gamma(m+d-1,m)=gf((d-1)/T+(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
% % % % %                 gamma(m,m+d-1)=subs(Gn,'f',(1-d)/T+(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
% % % % %                 gamma(m+d-1,m)=subs(Gn,'f',(d-1)/T+(cfo-param(i,2)*(13000+(m-1-K/2)/T))/(1+param(i,2)));
% % % %             end
% % % %         end
% % % %     end
% % % %     H=H+xind(i)*diag(exp(-sqrt(-1)*2*pi*param(i,3)/T*[-K/2:K/2-1]))*gamma;
% % % % end

% figure();
% image(abs(H),'CDataMapping','scaled');
% colorbar;

return

% % function x = OMP(A,b,sparsity)
% % %Step 1
% % index = []; k = 1; [Am, An] = size(A); r = b; x=zeros(An,1);
% % cor = A'*r; 
% % while k <= sparsity
% %     %Step 2
% %     [Rm,ind] = max(abs(cor)); 
% %     index = [index ind]; 
% %     %Step 3
% %     P = A(:,index)*inv(A(:,index)'*A(:,index))*A(:,index)';
% %     r = (eye(Am)-P)*b; 
% %     cor=A'*r;
% %     k=k+1;
% % end
% % %Step 5
% % xind = inv(A(:,index)'*A(:,index))*A(:,index)'*b;
% % x(index) = xind;
% % end

% % function theta=OMP(y,A,t)
% % [y_rows,y_columns] = size(y);
% % [M,N]=size(A);
% % theta=zeros(N,1);
% % At=[];
% % Pos_theta=[];
% % r_n=y;
% % product=A'*r_n;
% % while product>0.3
% %     [val,pos]=max(abs(product));
% %     At=[At,A(:,pos)];
% %     Pos_theta=[Pos_theta,pos];
% %     A(:,pos)=zeros(M,1);
% %     theta_ls=inv(At(:,1:end)'*At(:,1:end))
% % end