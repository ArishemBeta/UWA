function S_est=MIMO_LMMSE_Equalization(z,H,K,sc,Nr,Nt,MMode,nKr,N0,D)

z=reshape(z,Nr*K,1);
S_est=(H'*H+N0*eye(Nt*K))\H'*z;
S_est=reshape(S_est,K,Nt).';

% z=z.';
% if(strcmp(MMode,'QPSK'))
%     S_bar=zeros(Nt,K)
%     S_var=repmat(sc,Nt,K/length(sc));
% elseif(strcmp(MMode,'16QAM'))
%     S_bar=zeros(Nt,K);
%     S_var=repmat(sc,Nt,K/length(sc))*0.9472;
% end
% 
% S_est=zeros(Nt,K);
% noise=N0;
% 
% for m=1:K
%     Hm=H(:,m:K:end);
% %     Hm(Hm==0)=[];
%     Zm=reshape(z(:,max(m-D,1):min(m+D,K)),[],1);
% %     reshape(z,1,K*Nr);
% 
%     f=Hm*diag(S_var(:,m))*Hm'+noise*eye(height(Hm));
%     f=f\eye(size(f));
%     for nt=1:Nt
%         Hm_nt=Hm(:,nt);
%         S_est(nt,m)=Hm_nt'*f*(Zm-Hm*S_bar(:,m));
%     end
% end


% if(strcmp(MMode,'QPSK'))
%     S_bar=gpuArray.zeros(Nt*K,1);
%     S_var=gpuArray(diag(repmat(sc,1,Nt*K/length(sc))));
% elseif(strcmp(MMode,'16QAM'))
%     S_bar=gpuArray.zeros(Nt*K,1);
%     S_var=gpuArray(diag(repmat(sc,1,Nt*K/length(sc))*0.9472));
% end
% 
% S_est=gpuArray.zeros(Nt,K);
% z=gpuArray(reshape(z,Nr*K,1));
% % noise=sum(abs((diag(repmat(not(sc),1,K*Nr/length(sc)))*z)).^2)/(K*nKr*Nr);
% noise=N0;
% H=gpuArray(H);
% f=H*S_var*H'+noise*eye(K*Nr);
% f=f\eye(size(f));
% 
% for nt=1:Nt
%     for m=1:K
%         f_m=f*H(:,(nt-1)*K+m);
%         S_est(nt,m)=S_bar((nt-1)*K+m)+S_var((nt-1)*K+m,(nt-1)*K+m)...
%             *f_m'...
%             *(z-H*S_bar);
%     end
% end
% S_est=gather(S_est);