%Created by Jun Tao
%
function PlotChannel_new(H,Nt,Nr,L)

Nt
Nr
L

figure(101);

n= 1; k= 1;
subplot(2,2,1);        
plot(1:L,abs(H(k,(n-1)*L+1:n*L)));
xlabel('L');
title(['CIR of Channel: ',int2str(n),int2str(k)]);
grid;

n= 1; k= 8;
subplot(2,2,2);        
plot(1:L,abs(H(k,(n-1)*L+1:n*L)));
xlabel('L');
title(['CIR of Channel: ',int2str(n),int2str(k)]);
grid;

n= 2; k= 1;
subplot(2,2,3);        
plot(1:L,abs(H(k,(n-1)*L+1:n*L)));
xlabel('L');
title(['CIR of Channel: ',int2str(n),int2str(k)]);
grid;

n= 2; k= 8;
subplot(2,2,4);        
plot(1:L,abs(H(k,(n-1)*L+1:n*L)));
xlabel('L');
title(['CIR of Channel: ',int2str(n),int2str(k)]);
grid;