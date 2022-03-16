function UWAchannel=UWAchannel_generation(Npath,T,dT,D,dD,delay,doppler)
% Npath=5;
% T=0.2;dT=0.00025;
% D=80;dD=0.1;
% NAmplitude=4;
NDelay=1;
A=0.01;
diffD=[delay(1),doppler(1);
    delay(2),doppler(2);
    delay(3),doppler(3)];

t=[0:dT:T-dT];
d=[0:dD:D-dD];
UWAchannel=zeros(length(t),2*Npath);
% UWAchannel_plot=zeros(length(t),length(d));
for np=1:Npath
    tempA=A*1*0.5^(np-1)+0.00001*rand();
    tempD=diffD(np,1);
%     UWAchannel_plot(1,round(tempD/dD)+1)=tempA;
    UWAchannel(1,2*np-1)=tempA;
    UWAchannel(1,2*np)=tempD/1000;
end
for nt=2:length(t)
    for np=1:Npath
        tempA=UWAchannel(nt-1,2*np-1)*1+0.00001*(rand()-0.5);
        tempD=diffD(np,1)-diffD(np,2)*(nt-1)*dT*1000;
%          UWAchannel(nt,:)=UWAchannel(nt,:)+tempA*1*sign(dirac(d-tempD));
%         UWAchannel_plot(nt,round(tempD/dD)+1)=tempA;
        UWAchannel(nt,2*np-1)=tempA;
        UWAchannel(nt,2*np)=tempD/1000;
    end
end
% image(d,t,UWAchannel_plot,'CDataMapping','scaled');%
% colorbar;
% colormap('jet');
% caxis('auto');%[0,ceil(max(max(UWAchannel_plot)))]
% xlabel('delay(ms)','FontSize',12);
% ylabel('time(s)','FontSize',12);
return