function S_est=MIMO_LMMSE_Equalization(z,H,K,sc,Nr,Nt,MMode,nKr,N0)

if(strcmp(MMode,'QPSK'))
    S_bar=gpuArray.zeros(Nt*K,1);
    S_var=gpuArray(diag(repmat(sc,1,Nt*K/length(sc))));
elseif(strcmp(MMode,'16QAM'))
    S_bar=gpuArray.zeros(Nt*K,1);
    S_var=gpuArray(diag(repmat(sc,1,Nt*K/length(sc))*0.9472));
end

S_est=gpuArray.zeros(Nt*K,1);
z=gpuArray(reshape(z,Nr*K,1));
% noise=sum(abs((diag(repmat(not(sc),1,K*Nr/length(sc)))*z)).^2)/(K*nKr*Nr);
noise=N0;
H=gpuArray(H);
f=H*S_var*H'+noise*eye(K*Nr);
f=f\eye(size(f));
for m=1:K
    f_m=f*H(:,m:K:end);
    S_est((m-1)*Nt+1:m*Nt)=S_bar(m:K:end)+S_var(m:K:end,m:K:end)...
                            *f_m'...
                            *(z-H*S_bar);
end
S_est=gather(S_est);