function x = OMP(A,b,sparsity)
%Step 1
index = []; k = 1; [Am, An] = size(A); r = b; x=zeros(An,1);
cor = A'*r; 
while k <= sparsity
    %Step 2
    [Rm,ind] = max(abs(cor)); 
    index = [index ind]; 
    %Step 3
    P = A(:,index)*inv(A(:,index)'*A(:,index))*A(:,index)';
    r = (eye(Am)-P)*b; cor=A'*r;
    k=k+1;
end
%Step 5
xind = inv(A(:,index)'*A(:,index))*A(:,index)'*b;
x(index) = xind;
end