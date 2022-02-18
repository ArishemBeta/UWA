function UWAchannel=UWAchannel_generation(Npath,T,dT,D,dD,NAmplitude,NDelay)
T=2;dT=0.02;D=80;dD=0.1;Npath=5;NAmplitude=4;NDelay=1;
diffA=[9,0,-1,0;
    6,0,1,0;
    2,0,1,0;
    1,0,1,0;
    1,0,1,0];
diffD=[5,5;
    10,5;
    20,5;
    40,5;
    55,5];

t=[0:dT:T-dT];
d=[0:dD:D-dD];
UWAchannel=zeros(length(t),length(d));
for nt=1:T/dT
    for np=1:Npath
        tempA=0;
        for na=1:NAmplitude
            tempA=tempA+1/factorial(na-1)*diffA(np,na)*((nt-1)*dT)^(na-1);
        end
%         tempD=0;
%         for nd=1:NDelay
%             tempD=tempD+(-1)^nd/fractorial(nd)*diffD(nd)*((nt-1)*dT)^nd;
%         end
        tempD=diffD(np,1)+diffD(np,2)*(nt-1)*dT;
%          UWAchannel(nt,:)=UWAchannel(nt,:)+tempA*1*sign(dirac(d-tempD));
        UWAchannel(nt,round(tempD/dD))=tempA*0.1;
    end
end
image(d,t,UWAchannel,'CDataMapping','scaled');%
colorbar;
colormap('jet');
%  caxis([0,ceil(max(max(UWAchannel)))]);
xlabel('delay(ms)','FontSize',12);
ylabel('time(s)','FontSize',12);
return