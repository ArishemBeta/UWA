function [doppler_scalee,shifte]=mseqcorrelate_d_s(dl,dh,sl,sh,Nt,Nr,Ns,K,L,Kg,RX_data,nblk,doppler_scale,seq)
d=[dl:0.0001:dh];
co=zeros(1,length(d));
Noffset=1300*Nt+511+189+(nblk-1)*(K+2*Kg+length(seq));
Noffset=Noffset-fix(2000*doppler_scale(1));
if(nblk>1)
    for i=2:nblk
        Noffset=Noffset-floor(doppler_scale(i-1)*(K+2*Kg+length(seq)));
    end
    RX_block=RX_data(:,Noffset*Ns+1:Ns:(Noffset+length(seq))*Ns);
else
    RX_block=RX_data(:,Noffset*Ns+1:Ns:(Noffset+length(seq))*Ns);
end
% RX_block=(RX_data(:,Noffset+1:Ns:Noffset+length(seq)*Ns));
for k=1:length(d)
    Rx_block=RX_block;
    [I,D]=numden(sym(rats(d(k))));
    I=eval(I);D=eval(D);
    for i=1:Nr
%         Rx_block(i,:)=interp1((0:length(RX_block(i,:))-1),RX_block(i,:),(0:length(RX_block(i,:))-1)/(1+d(k)),'spline');
        seq1=resample(seq,D,D+I);
        for n=1:length(Rx_block(i,:))
            Rx_block(i,n)=Rx_block(i,n)*(exp(sqrt(-1)*2*pi*d(k)*(-13000)*(n-1)/9765.6));%1748828.125
        end
%         co(k,i,:)=xcorr(Rx_block(i,:),conj(seq));
%         co(k,(i-1)*(2*length(seq)-1)+1:i*(2*length(seq)-1))=xcorr(Rx_block(i,:),conj(seq));
%         co(k,:)=co(k,:)+abs(xcorr(Rx_block(i,:),conj(seq)));
%         co(k)=co(k)+sum(abs(Rx_block(i,:).*conj(seq)));
    end
    
%     co(k)=abs(co(k));
    co1=abs(xcorr(Rx_block(1,:),conj(seq1)));
    [value(k),index(k)]=max(co1);
    plot(co1);
    hold on;
    clear seq1;
end
plot(d,value);
hold on;
[val,idx]=max(value);
doppler_scalee=d(idx);
shifte=0;
return