function cfo=cfo_nullsubcarrier(range,Rx_block,K,sc,B,Ns)
ind_null=diag(not(repmat(sc,1,K/length(sc))));
[Nr, Len]=size(Rx_block);
e=zeros(1,length(range));
Rx_blockt=Rx_block;
for i=1:length(range)
    for j=1:Nr
        for n=1:Len
            Rx_blockt(j,n)=Rx_block(j,n)*(exp(sqrt(-1)*2*pi*(-range(i))*(n-1)/(B*Ns)));
        end
    end
    z=fft(Rx_blockt(:,1:Ns:K*Ns).');
    e(i)=sum(sum(abs((ind_null*z).^2)));
end
[me, ind_me]=min(e);
cfo=range(ind_me);
end