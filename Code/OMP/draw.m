Nr=2:2:12;
LS_LLR_Tx1=[0.3369	0.127	0.0674	0.0605	0.0557	0.0459];
LS_LLR_Tx2=[0.3477	0.1084	0.042	0.0713	0.0625	0.0498];
LS_LLR=(LS_LLR_Tx1+LS_LLR_Tx2)/2;
LS_Tx1=[0.376	0.2305	0.1377	0.1523	0.124	0.1328];
LS_Tx2=[0.3525	0.2393	0.1396	0.1465	0.1045	0.1426];
LS=(LS_Tx1+LS_Tx2)/2;
OMP_Tx1=[0.0371	0.0088	0	0	0	0];
OMP_Tx2=[0.0449	0	0	0	0	0];
OMP=(OMP_Tx1+OMP_Tx2)/2;

plot(Nr,LS_LLR);
hold on;
plot(Nr,LS);
plot(Nr,OMP);
xlabel('Phones','FontSize',20);
ylabel('BER','FontSize',20);
legend('LS\_LLR','LS','OMP');

plot(Nr,LS_LLR_Tx1);
hold on;
plot(LS_LLR_Tx2);
plot(LS_Tx1);
plot(LS_Tx2);
plot(OMP_Tx1);
plot(OMP_Tx2);