function UWAchannel=UWAchannel_generation(Npath,T,dT,D,dD)
% Npath=5;
% T=0.2;dT=0.00025;
% D=80;dD=0.1;
% NAmplitude=4;
NDelay=1;
A=0.01;
% diffA=[9,0,-1,0;
%     6,0,1,0;
%     2,0,1,0;
%     1,0,1,0;
%     1,0,1,0];
diffD=[2,24;
    8,24;
    59,30];

t=[0:dT:T-dT];
d=[0:dD:D-dD];
UWAchannel=zeros(length(t),2*Npath);
UWAchannel_plot=zeros(length(t),length(d));
for np=1:Npath
    tempA=A*1*0.5^(np-1)+0.00001*rand();
    tempD=diffD(np,1);
    UWAchannel_plot(1,round(tempD/dD))=tempA;
    UWAchannel(1,2*np-1)=tempA;
    UWAchannel(1,2*np)=tempD;
end
for nt=2:length(t)
    for np=1:Npath
%         tempA=0;
%         for na=1:NAmplitude
%             tempA=tempA+1/factorial(na-1)*diffA(np,na)*((nt-1)*dT)^(na-1);
            
%         end
%         tempD=0;
%         for nd=1:NDelay
%             tempD=tempD+(-1)^nd/fractorial(nd)*diffD(nd)*((nt-1)*dT)^nd;
%         end
        tempA=UWAchannel(nt-1,2*np-1)*1+0.00001*(rand()-0.5);
        tempD=diffD(np,1)+diffD(np,2)*(nt-1)*dT;
%          UWAchannel(nt,:)=UWAchannel(nt,:)+tempA*1*sign(dirac(d-tempD));
        UWAchannel_plot(nt,round(tempD/dD))=tempA;
        UWAchannel(nt,2*np-1)=tempA;
        UWAchannel(nt,2*np)=tempD;
    end
end
image(d,t,UWAchannel_plot,'CDataMapping','scaled');%
colorbar;
colormap('jet');
caxis('auto');%[0,ceil(max(max(UWAchannel_plot)))]
xlabel('delay(ms)','FontSize',12);
ylabel('time(s)','FontSize',12);
return