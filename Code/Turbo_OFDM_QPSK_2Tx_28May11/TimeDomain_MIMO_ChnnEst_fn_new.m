%-----------------------------------------------------
%  Time-Domain Channel Estimation
% 
%-----------------------------------------------------

function h = TimeDomain_MIMO_ChnnEst_fn_new(y, s, Nt, Nr, L, Np, SNR)
%it is noted that SNR here is refer to 'sigma_z^2/sigma_h^2' with
%'sigma_z^2' being the variance of noise and 'sigma_h^2'being the variance
%of single channel tap h_{m,n}(l)

h=zeros(Nr,L*Nt);
X=[];
for n=1: Nt
    if(0)
        c=s(n,L:Np).';
        r=s(n,L:-1:1);
    elseif(0)
        c=s(n,1: Np).';
        r= [s(n,1),zeros(1,L-1)];
    else
        c= [s(n,1: Np).';zeros(L-1,1)];
        r= [s(n,1),zeros(1,L-1)];        
    end
    X=[X toeplitz(c,r)];
end
    
for m=1:Nr
    if(0)
        Y=y(m, L:Np).';
    elseif(0)
        Y=y(m, 1:Np).';
    else
        Y=y(m, 1:Np+L-1).';
    end
    if(0)
        htmp=pinv(X'*X + 1/SNR*eye(L*Nt)) * X'* Y;
    elseif(0)       
        htmp=pinv(X'*X + 100*Nt/SNR*eye(L*Nt)) * X'* Y;
    else
        htmp=pinv(X'*X) * X'* Y;
    end
    h(m,:)=htmp.';
end

return
