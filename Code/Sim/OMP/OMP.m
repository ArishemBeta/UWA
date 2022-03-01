function H=OMP(z,Npa,Namp,B,Ns,K,sc_idx,block_symbol)
sparsity=Npa*((Namp+1)+1+1);
pilot_symbol=block_symbol(sc_idx).';
delay_hat=[0:1/(B*Ns):240/B];
beta=[-bmax:dbeta:bmax];
Np=length(delay_hat);                             %----------Nb-----------
Nb=length(beta);                                  %|
A=zeros(K/length(sc_idx),(Namp+1)*Np*Nb);         %|
Ablock=zeros(K/length(sc_idx),Np*Nb);             %Np
for n=1:Namp+1                                    %|
    for p=1:Np
        for b=1:Nb
            A(:,(n-1)*Np*Nb+(p-1)*Nb+b)=diag(e^(-sqrt(-1)*2*pi*delay_hat(p))*B/K*[-K/2:K/2-1])...
                ...
                *pilot_symbol;
        end
    end
end
xi=zeros((Namp+1)*Np*Nb,1);
% for n=1:Namp+1
%     for p=1:Np
%         xi((n-1)*Np+p)=(-1)^(n-1)/factorial(n-1)*(sqrt(-1)/(2*pi))^(n-1)*beta()
%     end
% end



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