%Compute the a priori statistics of symbols based on bit a priori LL
%Created by Jun Tao at MST on 24June11
%
function [S_bar S_var]= PrioriStatCal_QPSK_fn(Nt,K,LLa_cod);

Q= 4; 
QPSK_Sym_Set= [1 -i i -1]; 

LLa_sym= zeros(Q*Nt,K);
norm_fac= zeros(1,K);
LLa_sym_norm= zeros(Q,K);
Pa_sym= zeros(Q,K);

S_bar= zeros(Nt, K); 
S_var= zeros(Nt, K); 

for nt= 1: Nt
    %convert bit LL to symbol LL
    LLa_sym((nt-1)*Q+1,:)= LLa_cod((nt-1)*2+1,1:2:end)+LLa_cod((nt-1)*2+1,2:2:end); %1
    LLa_sym((nt-1)*Q+2,:)= LLa_cod((nt-1)*2+2,1:2:end)+LLa_cod((nt-1)*2+1,2:2:end); %-i
    LLa_sym((nt-1)*Q+3,:)= LLa_cod((nt-1)*2+1,1:2:end)+LLa_cod((nt-1)*2+2,2:2:end); %i
    LLa_sym((nt-1)*Q+4,:)= LLa_cod((nt-1)*2+2,1:2:end)+LLa_cod((nt-1)*2+2,2:2:end); %-1
    %normalize symbol a priori probability
    for m= 1: K
        norm_fac(m)= -logsum(LLa_sym((nt-1)*Q+1: nt*Q,m).');
    end
    LLa_sym_norm= LLa_sym((nt-1)*Q+1: nt*Q,:)+ ones(Q,1)*norm_fac;
    Pa_sym= exp(LLa_sym_norm);
    Pa_sym= Pa_sym./(ones(Q,1)*sum(Pa_sym,1));
    
    %calculate the mean and variance of the symbols   
    S_bar(nt,:)= QPSK_Sym_Set* Pa_sym; 
    %calculate a priori variance of Tx symbols
    S_var(nt,:)= abs(QPSK_Sym_Set).^2* Pa_sym- abs(S_bar(nt,:)).^2;
end