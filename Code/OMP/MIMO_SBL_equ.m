function S_est = MIMO_SBL_equ(z,H,r_max,error,noise,sc,Nt,K,Var)
%z:Nr*Kά�ӽ���������
%H:(Nr*K,Ht*K)ά�ŵ�����
%r_max:����������(�Ƽ����ã�10)
%error:�жϵ�������������ֵ���Ƽ����ã�1e-6)
%noise:��������
%S_var:�������ŷ���ľ���
z=reshape(z,[],1);
S_var=zeros(1,Nt*K);
for nt=1:Nt
    tv=zeros(1,K);
    t=repmat(sc(nt,:),1,K/width(sc));
    tv=tv+(t==2)*Var(nt,2)+(t==1)*Var(nt,1);
    S_var((nt-1)*K+1:nt*K)=tv;
end
sigma_S = diag(S_var);
for i = 0:r_max
% E-step
cov = inv((1/noise)*H'*H+S_var);
mu = (1/noise)*cov*H'*z;
%������һ�ε����õ��ķ��ŷ���
sigma_S_store = sigma_S;
%M-step
    sigma_S = diag(cov)+abs(mu).^2;
a =  norm(sigma_S-sigma_S_store,2)^2;
if  norm(sigma_S-sigma_S_store,2)^2 < error
    break;
end
end 
S_est = reshape(mu,[],2).';