%Compute the a priori statistics of symbols based on bit a priori LL
%Created by Jun Tao at MST on 25June11
%
function [S_bar S_var]= PrioriStatCal_16QAM_fn(Nt,K,LLa_cod)

Q= 16; 
HQAM_Sym_Set= [-0.9487 - 0.9487i, -0.9487 - 0.3162i, -0.9487 + 0.9487i, -0.9487 + 0.3162i,...
               -0.3162 - 0.9487i, -0.3162 - 0.3162i, -0.3162 + 0.9487i, -0.3162 + 0.3162i,...
                0.9487 - 0.9487i,  0.9487 - 0.3162i,  0.9487 + 0.9487i,  0.9487 + 0.3162i,...
                0.3162 - 0.9487i,  0.3162 - 0.3162i,  0.3162 + 0.9487i,  0.3162 + 0.3162i]; 

LLa_sym= zeros(Q*Nt,K);
norm_fac= zeros(1,K);
LLa_sym_norm= zeros(Q,K);
Pa_sym= zeros(Q,K);

S_bar= zeros(Nt, K); 
S_var= zeros(Nt, K); 

for nt= 1: Nt
    %convert bit LL to symbol LL
    LLa_sym((nt-1)*Q+1,:)= LLa_cod((nt-1)*2+1,1:4:end)+LLa_cod((nt-1)*2+1,2:4:end)+LLa_cod((nt-1)*2+1,3:4:end)+LLa_cod((nt-1)*2+1,4:4:end); %sym 1
    LLa_sym((nt-1)*Q+2,:)= LLa_cod((nt-1)*2+1,1:4:end)+LLa_cod((nt-1)*2+1,2:4:end)+LLa_cod((nt-1)*2+1,3:4:end)+LLa_cod((nt-1)*2+2,4:4:end); %sym 2
    LLa_sym((nt-1)*Q+3,:)= LLa_cod((nt-1)*2+1,1:4:end)+LLa_cod((nt-1)*2+1,2:4:end)+LLa_cod((nt-1)*2+2,3:4:end)+LLa_cod((nt-1)*2+1,4:4:end); %sym 3
    LLa_sym((nt-1)*Q+4,:)= LLa_cod((nt-1)*2+1,1:4:end)+LLa_cod((nt-1)*2+1,2:4:end)+LLa_cod((nt-1)*2+2,3:4:end)+LLa_cod((nt-1)*2+2,4:4:end); %sym 4
    LLa_sym((nt-1)*Q+5,:)= LLa_cod((nt-1)*2+1,1:4:end)+LLa_cod((nt-1)*2+2,2:4:end)+LLa_cod((nt-1)*2+1,3:4:end)+LLa_cod((nt-1)*2+1,4:4:end); %sym 5
    LLa_sym((nt-1)*Q+6,:)= LLa_cod((nt-1)*2+1,1:4:end)+LLa_cod((nt-1)*2+2,2:4:end)+LLa_cod((nt-1)*2+1,3:4:end)+LLa_cod((nt-1)*2+2,4:4:end); %sym 6
    LLa_sym((nt-1)*Q+7,:)= LLa_cod((nt-1)*2+1,1:4:end)+LLa_cod((nt-1)*2+2,2:4:end)+LLa_cod((nt-1)*2+2,3:4:end)+LLa_cod((nt-1)*2+1,4:4:end); %sym 7
    LLa_sym((nt-1)*Q+8,:)= LLa_cod((nt-1)*2+1,1:4:end)+LLa_cod((nt-1)*2+2,2:4:end)+LLa_cod((nt-1)*2+2,3:4:end)+LLa_cod((nt-1)*2+2,4:4:end); %sym 8

    LLa_sym((nt-1)*Q+9,:)= LLa_cod((nt-1)*2+2,1:4:end)+LLa_cod((nt-1)*2+1,2:4:end)+LLa_cod((nt-1)*2+1,3:4:end)+LLa_cod((nt-1)*2+1,4:4:end); %sym 9
    LLa_sym((nt-1)*Q+10,:)= LLa_cod((nt-1)*2+2,1:4:end)+LLa_cod((nt-1)*2+1,2:4:end)+LLa_cod((nt-1)*2+1,3:4:end)+LLa_cod((nt-1)*2+2,4:4:end); %sym 10
    LLa_sym((nt-1)*Q+11,:)= LLa_cod((nt-1)*2+2,1:4:end)+LLa_cod((nt-1)*2+1,2:4:end)+LLa_cod((nt-1)*2+2,3:4:end)+LLa_cod((nt-1)*2+1,4:4:end); %sym 11
    LLa_sym((nt-1)*Q+12,:)= LLa_cod((nt-1)*2+2,1:4:end)+LLa_cod((nt-1)*2+1,2:4:end)+LLa_cod((nt-1)*2+2,3:4:end)+LLa_cod((nt-1)*2+2,4:4:end); %sym 12
    LLa_sym((nt-1)*Q+13,:)= LLa_cod((nt-1)*2+2,1:4:end)+LLa_cod((nt-1)*2+2,2:4:end)+LLa_cod((nt-1)*2+1,3:4:end)+LLa_cod((nt-1)*2+1,4:4:end); %sym 13
    LLa_sym((nt-1)*Q+14,:)= LLa_cod((nt-1)*2+2,1:4:end)+LLa_cod((nt-1)*2+2,2:4:end)+LLa_cod((nt-1)*2+1,3:4:end)+LLa_cod((nt-1)*2+2,4:4:end); %sym 14
    LLa_sym((nt-1)*Q+15,:)= LLa_cod((nt-1)*2+2,1:4:end)+LLa_cod((nt-1)*2+2,2:4:end)+LLa_cod((nt-1)*2+2,3:4:end)+LLa_cod((nt-1)*2+1,4:4:end); %sym 15
    LLa_sym((nt-1)*Q+16,:)= LLa_cod((nt-1)*2+2,1:4:end)+LLa_cod((nt-1)*2+2,2:4:end)+LLa_cod((nt-1)*2+2,3:4:end)+LLa_cod((nt-1)*2+2,4:4:end); %sym 16
       
    %normalize symbol a priori probability
    for m= 1: K
        norm_fac(m)= -logsum(LLa_sym((nt-1)*Q+1: nt*Q,m).');
    end
    LLa_sym_norm= LLa_sym((nt-1)*Q+1: nt*Q,:)+ ones(Q,1)*norm_fac;
    Pa_sym= exp(LLa_sym_norm);
    Pa_sym= Pa_sym./(ones(Q,1)*sum(Pa_sym,1));
    
    %calculate the mean and variance of the symbols   
    S_bar(nt,:)= HQAM_Sym_Set* Pa_sym; 
    %calculate a priori variance of Tx symbols
    S_var(nt,:)= abs(HQAM_Sym_Set).^2* Pa_sym- abs(S_bar(nt,:)).^2;
end

