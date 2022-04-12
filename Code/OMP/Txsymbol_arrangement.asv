function [Tx_sym_insert, Tx_sym]=Txsymbol_arrangement(data_sym,pilot_sym,K,pattern,Ns,Nt)

l=width(pattern);
g=K/l;
Tx_sym_insert=zeros(K*Ns,width(data_sym));
Tx_sym=zeros(K,width(data_sym));
d=1;
p=1;

for i=1:g
    for j=1:l
        if(pattern(j)==1)
            if((i-1)*l+j<K/2)
                Tx_sym_insert((i-1)*l+j,:)=data_sym(d,:);
            else
                Tx_sym_insert((i-1)*l+j+(Ns-1)*K,:)=data_sym(d,:);
            end
            Tx_sym((i-1)*l+j,:)=data_sym(d,:);
            d=d+1;
        elseif(pattern(j)==2)
            if((i-1)*l+j<K/2)
                Tx_sym_insert((i-1)*l+j,:)=pilot_sym(p,:);
            else
                Tx_sym_insert((i-1)*l+j+(Ns-1)*K,:)=pilot_sym(p,:);
            end
            Tx_sym((i-1)*l+j,:)=data_sym(p,:);
            p=p+1;
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