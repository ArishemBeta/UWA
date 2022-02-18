%Created by Jun Tao
%
function PlotChannel(H,Nt,N_chnn,L);

% Nt= 2;
% N_chnn= 12;
%L=100;
for n=1: Nt
%     figure(100+n*3);
%     clf;
    for k=1: N_chnn
        figure();      
        plot(1:L,abs(H(k,(n-1)*L+1:n*L)));
        xlabel('抽头数',"FontName","宋体","FontSize",16);
        ylabel('幅度',"FontName","宋体","FontSize",16);
        title(['冲激响应: ','接收机',int2str(k)],"FontName","宋体","FontSize",16);
        grid;
    end    
    
%     figure(100+n*3+1);
%     clf;
%     for k=N_chnn/2+1: N_chnn
%         subplot(2,2,k-N_chnn/2);        
%         plot(1:L,abs(H(k,(n-1)*L+1:n*L)));
%         xlabel('Tap index');
%         ylabel('Magnitude');
%         title(['CIR of Subchannel: ','T',int2str(n),' to H',int2str(k)]);
%         grid;
%     end     
end

%energy analysis for CIR
% for n=1: Nt
%     for k=1: N_chnn
%         E(n,k)= sum(abs(H(k,(n-1)*L+1: n*L)).^2);
%     end
% end
% E
% mean(E)*2
% Diff_E= E(1,:)-E(2,:)
% Diff_E_Ratio= E(1,:)./E(2,:)

% figure(1);
% plot(1: N_chnn, E(1,:), '-ro', 1: N_chnn, E(2,:), '-b*',1: N_chnn, Diff_E, '-ks');
% xlabel('Receiver Index');
% ylabel('CIR energy');
% legend('Tx1','Tx2','Tx1-Tx2');
% grid;