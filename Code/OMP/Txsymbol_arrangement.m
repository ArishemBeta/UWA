function [Tx_sym_insert, Tx_sym]=Txsymbol_arrangement(data_sym,pilot_sym,K,pattern,Ns,Nblock)

l=width(pattern);
g=K/l;
Nt=height(pattern);
Nblock=width(data_sym)/Nt;
Tx_sym_insert=zeros(K*Ns,Nblock*Nt);
Tx_sym=zeros(K,Nblock*Nt);

for nt=1:Nt
    d=1;
    p=1;
    for i=1:g
        for j=1:l
            if(pattern(nt,j)==1)
                if((i-1)*l+j<K/2)
                    Tx_sym_insert((i-1)*l+j,(nt-1)*Nblock+1:nt*Nblock)=data_sym(d,(nt-1)*Nblock+1:nt*Nblock);
                else
                    Tx_sym_insert((i-1)*l+j+(Ns-1)*K,(nt-1)*Nblock+1:nt*Nblock)=data_sym(d,(nt-1)*Nblock+1:nt*Nblock);
                end
                Tx_sym((i-1)*l+j,(nt-1)*Nblock+1:nt*Nblock)=data_sym(d,(nt-1)*Nblock+1:nt*Nblock);
                d=d+1;
            elseif(pattern(nt,j)==2)
                if((i-1)*l+j<K/2)
                    Tx_sym_insert((i-1)*l+j,(nt-1)*Nblock+1:nt*Nblock)=pilot_sym(p,(nt-1)*Nblock+1:nt*Nblock);
                else
                    Tx_sym_insert((i-1)*l+j+(Ns-1)*K,(nt-1)*Nblock+1:nt*Nblock)=pilot_sym(p,(nt-1)*Nblock+1:nt*Nblock);
                end
                Tx_sym((i-1)*l+j,(nt-1)*Nblock+1:nt*Nblock)=pilot_sym(p,(nt-1)*Nblock+1:nt*Nblock);
                p=p+1;
            end
        end
    end
end

% l=length(pattern);
% g=K/l;
% Tx_sym_insert=zeros(K*Ns,width(data_sym));
% Tx_sym=zeros(K,width(data_sym));
% p=1;
% 
% for i=1:g
%     for j=1:l
%         if(pattern(j)==1)
%             if((i-1)*l+j<K/2)
%                 Tx_sym_insert((i-1)*l+j,:)=data_sym(p,:);
%             else
%                 Tx_sym_insert((i-1)*l+j+(Ns-1)*K,:)=data_sym(p,:);
%             end
%             Tx_sym((i-1)*l+j,:)=data_sym(p,:);
%             p=p+1;
%         end
%     end
% end

return