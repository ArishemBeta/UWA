function [doppler_scalee,shifte]=mseqcorrelate_d_s(dl,dh,sl,sh,Nt,Nr,Ns,K,Kg,RX_data,nblk,doppler_scale,seq)
d=[dl:0.0001:dh];
co=zeros(1,length(d));
Noffset=(nblk-1)*(K+2*Kg+length(seq));
% doppler_scale=[0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005];
shift=0;%(nblk-1)*1
if(nblk>1)
    for i=2:nblk
        Noffset=Noffset-round(doppler_scale(i-1)*(K+2*Kg+length(seq)));
    end
    RX_block=RX_data(:,(Noffset+shift)*Ns+1:(shift+Noffset+round((length(seq))*(1-doppler_scale(nblk-1))))*Ns);
else
    RX_block=RX_data(:,Noffset*Ns+1:(Noffset+length(seq))*Ns);
end
% RX_block=RX_data(:,Noffset*Ns+1:(Noffset+length(seq))*Ns);
% for i=1:Nr
%     meanseq=mean(abs(RX_block(i,:)));
%     for j=1:length(RX_block(i,:))
%         if(abs(RX_block(i,j))<meanseq) RX_block(i,j)=0;
%         else RX_block(i,j)=0.015;
%         end
%     end
% end

% plot(abs(RX_block(1,:)),'.');
for k=1:length(d)
    [I,D]=numden(sym(rats(d(k))));
    I=eval(I);D=eval(D);
    clear Rx_block;
    co1=0;
    for i=1:Nr
%         Rx_block=RX_block;
%         Rx_block(i,:)=interp1((0:length(RX_block(i,:))-1),RX_block(i,:),(0:length(RX_block(i,:))-1)/(1+d(k)),'spline');
        Rx_block(i,:)=resample(RX_block(i,:),D+I,D);
%         seq1=resample(seq,D,D+I);
%         for n=1:length(Rx_block(i,:))
%             Rx_block(i,n)=Rx_block(i,n)*(exp(sqrt(-1)*2*pi*(-d(k))*4882*(n-1)*^2/9765.625));%48828.125
%         end(nblk-1)*(K+2*Kg+length(seq))+
%         for n=1:length(Rx_block(i,:))
%             Rx_block(i,n)=Rx_block(i,n)*(exp(sqrt(-1)*2*pi*(-d(k))*(n-1)*(I/(D+I))/(9765.625/2)));%48828.125
%         end
%         for n=1:length(Rx_block(i,:))
%             Rx_block(i,n)=Rx_block(i,n)*(exp(sqrt(-1)*2*pi*(d(k))*13000*(n-1)/9765.625));%48828.125
%         end
%         co(k)=co(k)+sum(abs(Rx_block(i,:).*conj(seq)));
        co1=co1+abs(xcorr(Rx_block(i,:),conj(seq)));
    end
    
    
    [value(k),index(k)]=max(co1);
%     value(k)=co1(127);
%     plot(co1);
%     hold on;
    
end
% plot(abs(Rx_block(1,:)),'.');
plot(d,value);
hold on;
[val,idx]=max(value);
doppler_scalee=d(idx);
% doppler_scalee=d(idx)+0.001*(nblk-1);
% if(doppler_scalee>0.005) doppler_scalee=0.005;end
shifte=0;
return

% function [doppler_scalee,shifte]=mseqcorrelate_d_s(dl,dh,sl,sh,Nt,Nr,Ns,K,L,Kg,RX_data,nblk,doppler_scale,seq)
% d=[dl:0.0001:dh];
% co=zeros(1,length(d));
% Noffset=1300*Nt+511+189+(nblk-1)*(K+2*Kg+length(seq));
% for i=1:nblk
%     Noffset=Noffset-floor(doppler_scale(max(1,i-1))*(K+2*Kg+length(seq)));
% end
% RX_block=(RX_data(:,Noffset+1:Ns:Noffset+length(seq)*Ns));
% for k=1:length(d)
%     for i=1:Nr
%         Rx_block(i,:)=interp1((0:length(RX_block(i,:))-1),RX_block(i,:),(0:length(RX_block(i,:))-1)/(1+d(k)),'spline');
%         for n=1:length(Rx_block(i,:))
%             Rx_block(i,n)=Rx_block(i,n)*(exp(sqrt(-1)*2*pi*d(k)*(13000)*(n-1)/9765.6));%1748828.125
%         end
% %         co(k,i,:)=xcorr(Rx_block(i,:),conj(seq));
% %         co(k,(i-1)*(2*length(seq)-1)+1:i*(2*length(seq)-1))=xcorr(Rx_block(i,:),conj(seq));
% %         co(k,:)=co(k,:)+abs(xcorr(Rx_block(i,:),conj(seq)));
%         co(k)=co(k)+sum(abs(Rx_block(i,:).*conj(seq)));
%     end
%     co(k)=abs(co(k));
% end
% plot(co);
% hold on;
% [value,idx]=max(co);
% % plot(co_max);
% doppler_scalee=d(idx);
% shifte=0;
% return