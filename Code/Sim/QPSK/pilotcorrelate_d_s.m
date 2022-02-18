function [doppler_scalee,shifte]=pilotcorrelate_d_s(dl,dh,sl,sh,sc_idx,Nt,Nr,Ns,K,L,RX_data,nblk,doppler_scale,block_symbol)
Kg=240;
d=[dl:0.0001:dh];
Noffset=(nblk-1)*(K+Kg);
for i=1:nblk
    Noffset=Noffset-floor(doppler_scale(max(1,i-1))*(K+Kg));
end
RX_block=RX_data(:,Noffset+1:Noffset+(K+L));
pilot_symbol=zeros(1,K);
pilot_symbol(sc_idx)=block_symbol(sc_idx);
pilot_t=ifft(pilot_symbol.').';

for k=1:length(d)
    Rx_block=RX_block;
    for i=1:Nr
        Rx_block(i,:)=interp1((0:length(Rx_block(i,:))-1),Rx_block(i,:),(0:length(Rx_block(i,:))-1)/(1+d(k)),'spline');
    end
    co(k,:)=xcorr(Rx_block(1,:),conj(pilot_t));
    plot(abs(co(k,:)));
    hold on;
    co_max(k)=max(abs(co(k,:)));
end
[value,idx]=max(co_max);
doppler_scalee=d(idx);
shifte=0;
return