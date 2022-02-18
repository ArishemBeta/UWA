
function SoftSym= SoftSymCal_QPSK(LLa_cod,Nt,Ndata)

QPSK_Sym_Set= [1 -i i -1]; 
QPSK_Bit_Set= [0 0;1 0;0 1;1 1];
Q= 4; %four symbols in constellation for QPSK

%-----soft-decision symbols after channel MAP decoder----
for nt= 1: Nt
    %convert bit LL to symbol LL
    LLa_sym((nt-1)*Q+1,:)= LLa_cod((nt-1)*2+1,1:2:end)+LLa_cod((nt-1)*2+1,2:2:end);
    LLa_sym((nt-1)*Q+2,:)= LLa_cod((nt-1)*2+2,1:2:end)+LLa_cod((nt-1)*2+1,2:2:end);
    LLa_sym((nt-1)*Q+3,:)= LLa_cod((nt-1)*2+1,1:2:end)+LLa_cod((nt-1)*2+2,2:2:end);
    LLa_sym((nt-1)*Q+4,:)= LLa_cod((nt-1)*2+2,1:2:end)+LLa_cod((nt-1)*2+2,2:2:end);

    %----Calculate APP of current symbol being equalized conditioned on Rx vector----       
    %normalize symbol a priori probability
    for m= 1: Ndata
        norm_fac(m)= -logsum(LLa_sym((nt-1)*Q+1: nt*Q,m).');
    end
    LLa_sym_norm= LLa_sym((nt-1)*Q+1: nt*Q,:)+ ones(Q,1)*norm_fac;

    SoftSym(nt,:)= QPSK_Sym_Set*exp(LLa_sym_norm);
end