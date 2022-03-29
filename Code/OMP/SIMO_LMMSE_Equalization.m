function S_est=SIMO_LMMSE_Equalization(z,H,K,sc,Nr,Nt,MMode,nKr,N0)
if(strcmp(MMode,'QPSK'))
    S_bar=zeros(K,1);
    S_var=diag(repmat(sc,1,K/length(sc)));
elseif(strcmp(MMode,'16QAM'))
    S_bar=zeros(K,1);
    S_var=diag(repmat(sc,1,K/length(sc))*0.9472);
end

S_est=zeros(K,1);
z=reshape(z,Nr*K,1);
% noise=sum(abs((diag(repmat(not(sc),1,K*Nr/length(sc)))*z)).^2)/(K*nKr*Nr);
noise=N0;
f=inv(H*S_var*H'+noise*eye(K*Nr));
for m=1:K
%     m
    f_m=f*H(:,m);
    S_est(m)=S_bar(m)+S_var(m,m)*f_m'*(z-H*S_bar);
end