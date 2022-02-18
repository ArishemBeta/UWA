%-----------------------------------------------------
%  Time-Domain Channel Estimation
% 
%  Modified by Chengshan Xiao on October 18, 2008
%-----------------------------------------------------

function h = TimeDomain_MIMO_ChnnEst_fn(y, s, Nt, Nr, L, Np, SNR)

h=zeros(Nr,L*Nt);
X=[];
for n=1: Nt
    c=s(n,L:Np).';
    r=s(n,L:-1:1);
    X=[X toeplitz(c,r)];
end
    
for m=1:Nr
    Y=y(m, L:Np).';
    htmp=pinv(X'*X + 1/SNR*eye(L*Nt)) * X'* Y;
    h(m,:)=htmp.';
end

return
