function [Tx_sym_insert, Tx_sym]=Txsymbol_arrangement(data_sym,K,sc,Ns)
l=length(sc);
g=K/l;
Tx_sym_insert=zeros(K*Ns,width(data_sym));
Tx_sym=zeros(K,width(data_sym));
p=1;

for i=1:g
    for j=1:l
        if(sc(j)==1)
            if((i-1)*l+j<K/2)
                Tx_sym_insert((i-1)*l+j,:)=data_sym(p,:);
            else
                Tx_sym_insert((i-1)*l+j+(Ns-1)*K,:)=data_sym(p,:);
            end
            Tx_sym((i-1)*l+j,:)=data_sym(p,:);
            p=p+1;
        end
    end
end

return