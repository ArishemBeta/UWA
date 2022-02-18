a=[-5,-4,-5,-5,-5,-6,-5,-6];
plot(a,'-*','LineWidth',1.5,'MarkerSize',10);
ylabel("同步值（单位：样本点）","FontName","宋体","FontSize",16);
xlabel("块编号","FontName","宋体","FontSize",16);

a=[-0.06,-0.05,0.04,-0.04,0.08,-0.08,-0.09,0.01];
plot(a,'-*','LineWidth',1.5,'MarkerSize',10);
ylabel("CFO (Hz)","FontName","宋体","FontSize",16);
xlabel("Block Index","FontName","宋体","FontSize",16);

a=[-1e-06,1.0e-05,-1.5e-05,-0.5e-05,-1.2e-05,-1.0e-05,1.0e-05,1.2e-05];
plot(a,'-*','LineWidth',1.5,'MarkerSize',10);
ylabel("多普勒缩放因子","FontName","宋体","FontSize",16);
xlabel("块编号","FontName","宋体","FontSize",16);

a=[0.5,0.501,0.4824,0.5127,0.5020,0.4404,0.4883,0.5029,0.4893,0.4883,0.4844,0.4512];
b=[0.4795,0.4727,0.5078,0.499,0.5029,0.4873,0.4697,0.4980,0.4961,0.4883,0.4756,0.4873];
c=[0,0,0,0,0,0,0,0,0,0,0,0];
plot(a,'-+','LineWidth',1.5,'Color','r','MarkerSize',10,'MarkerEdgeColor','r');hold on;
plot(b,'-o','LineWidth',1.5,'Color','g','MarkerSize',10,'MarkerEdgeColor','g');hold on;
plot(c,'-*','LineWidth',1.5,'Color','b','MarkerSize',10,'MarkerEdgeColor','b');hold on;
ylabel("误码率","FontName","宋体","FontSize",16);
xlabel("块编号","FontName","宋体","FontSize",16);

figure();
subplot(3,1,1);
a=[10,10,10,10,20,20,20,20,10,10,20,20];
plot(a,'-*','LineWidth',1.5,'MarkerSize',10);hold on;
ylabel('$\hat{\tau}$ $(1/f_{s})$','interpreter','latex',"FontName","Times New Roman","FontSize",30);
xlabel("Block Index","FontName","Times New Roman","FontSize",24);
subplot(3,1,2);
b=[-0.0013985,-0.001499,-0.001298,-0.001499,-0.001298,-0.001499,-0.001298,-0.001097,-0.001298,-0.001499,-0.0017,-0.001097];
plot(b,'-*','LineWidth',1.5,'MarkerSize',10);hold on;
ylabel('$\hat{a}$','interpreter','latex',"FontName","Times New Roman","FontSize",30);
xlabel("Block Index","FontName","Times New Roman","FontSize",24);
subplot(3,1,3);
c=[-19.003,-19.003,-19.003,-18.7015,-19.003,-19.003,-19.3045,-19.3045,-19.003,-19.003,-18.4,-19.3045];
plot(c,'-*','LineWidth',1.5,'MarkerSize',10);hold on;
ylabel('$\hat{\varepsilon}$ (Hz)','interpreter','latex',"FontName","Times New Roman","FontSize",30);
xlabel("Block Index","FontName","Times New Roman","FontSize",24);

set(gca,'XTick',[1:1:12]);
set(gca,'YTick',[17.5:0.5:20.5]);

a=[0,0,0,0,0,0,0,0,0,0,0,0];%3
b=[0.5,0.501,0.4824,0.5127,0.5020,0.4404,0.4883,0.5029,0.4893,0.4883,0.4844,0.4512];%
c=[0.4727,0.4883,0.4688,0.4883,0.4883,0.4951,0.4717,0.5010,0.5225,0.4785,0.5264,0.4795];%at
d=[0.2109,0.2002,0.2305,0.2617,0.2617,0.2344,0.2754,0.2324,0.2529,0.2324,0.1738,0.2402];%ct
e=[0,0.001,0,0,0.0029,0,0,0,0,0,0.001,0.002];%ac
plot(a,'-*','LineWidth',1.5,'MarkerSize',10);hold on;
plot(b,'-o','LineWidth',1.5,'MarkerSize',10);hold on;
plot(c,'-+','LineWidth',1.5,'MarkerSize',10);hold on;
plot(d,'--','LineWidth',1.5,'MarkerSize',10);hold on;
plot(e,'-.','LineWidth',1.5,'MarkerSize',10);hold on;
ylabel("误码率","FontName","宋体","FontSize",20);
xlabel("块编号","FontName","宋体","FontSize",20);
h=legend('Scheme 1:$\tau,a,\varepsilon$','Scheme 2:$null$','Scheme 3:$\tau,a$','Scheme 4:$\tau,\varepsilon$','Scheme 5:$a,\varepsilon$');
set(h,'Interpreter','latex');

a=[0.030,0.025,0.042,0.038,0.055,0.096,0.054,0.041];%3
b=[0.0674,0.0430,0.1016,0.0752,0.1152,0.1924,0.1211,0.0830];%a
c=[0.0605,0.0498,0.0869,0.0820,0.1162,0.1924,0.1104,0.0820];%c
%d=[0.0674,0.0508,0.0986,0.0820,0.1230,0.1973,0.1172,0.0830];%t
e=[0.032,0.028,0.042,0.045,0.061,0.102,0.054,0.045];%ac
plot(a,'-*','LineWidth',1.5,'MarkerSize',10);hold on;
plot(b,'-o','LineWidth',1.5,'MarkerSize',10);hold on;
plot(c,'-+','LineWidth',1.5,'MarkerSize',10);hold on;
%plot(d,'--','LineWidth',1.5,'MarkerSize',10);hold on;
plot(e,'-.','LineWidth',1.5,'MarkerSize',10);hold on;
ylabel("误码率","FontName","宋体","FontSize",16);
xlabel("块编号","FontName","宋体","FontSize",16);
h=legend('Scheme 1:$\tau,a,\varepsilon$','Scheme 3:$\tau,a$','Scheme 4:$\tau,\varepsilon$','Scheme 5:$a,\varepsilon$');
set(h,'Interpreter','latex');

a=[1,1,1,1,1,1];
b=[1,1.2,3,3.6,5,6];
c=[0,0,0,0,0,0];
plot(b,a,'^');

a=[1,1,1];
b=[1,3,5,];
c=[0,0,0,0,0,0];
plot(b,a,'^');

a=[1,1,1,1,1,1];
b=[1.1,1.2,3.5,3.6,5.9,6];
c=[0,0,0,0,0,0];
plot(b,a,'^');

%[0.005859375,0.0029296875,0,0,0,0,0,0,0.001953125,0,0,0]
a=[0.05, 0.00138, 0.00065, 0];
b=[2,4,6,8];
plot(b,a,'-*','LineWidth',1.5,'MarkerSize',10);hold on;
ylabel("BER","FontName","Times New Roman","FontSize",24);
xlabel("Number of hydrophone","FontName","Times New Roman","FontSize",24);
% h=legend('Scheme 1:$\tau,a,\varepsilon$','Scheme 2:$null$','Scheme 3:$\tau,a$','Scheme 4:$\tau,\varepsilon$','Scheme 5:$a,\varepsilon$');
% set(h,'Interpreter','latex');
%8 0
%7
%6 0.00065
%5
%4 0.00138
%3
%2 0.05

