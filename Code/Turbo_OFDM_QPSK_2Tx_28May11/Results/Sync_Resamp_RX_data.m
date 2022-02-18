% ---Created by J.Tao at UMC on Dec.23, 2007 for 2006 moving data---
%
function RX_resample_data= Sync_Resamp_RX_data(RX_baseband,Ngap,Noff,Nr,Nt,Nsps1,N_chnn);

% --------- decimate the signal to have 2 samples per symbol-----------
Nsps=2;        %Number of samples per symbol %must be 1, 2, or 4
%N_chnn= 8;
Ndata=(11.9*4000)*Nsps;  %=47600*Nsps; in symbols, J.Tao
RX_data= zeros(Ndata,N_chnn);

if(1)
    Ngap=Ngap+Noff;  %introduce an offset to adjust current synchronization, J.Tao, Nov.17, 07
end
if(0)
    for n=1:N_chnn
        for k=1:Ndata
            RX_data(k,n)=RX_baseband(Ngap(n)+(k-1)*(20/Nsps),n); 
        end
    end
else
    for n=1:N_chnn
        RX_data(:,n)=RX_baseband(Ngap(n)+((1:Ndata)-1)*(20/Nsps),n); 
    end
end
% --------- resample Rx data-----------
P=95200/(20/Nsps)/2; %/4;
Q=round(Nr/Nt*P); %Q is a vector since Nr is a vector, J. Tao
%Q=round(mean(Q));

RX_resample_data=zeros(max(ceil(Ndata*P./Q)),N_chnn);    
for Chnn=1:N_chnn
    temp= resample(RX_data(:,Chnn), P, Q(Chnn));
    RX_resample_data(1:length(temp),Chnn)= temp;
end
if Nsps1==1
    if(1) %obtain one-sample-per-symbol sequence 
        %by downsampling the generated two-samples-per-symbol sequence
        RX_resample_data= RX_resample_data(1:2:end,:); 
    else %obtain one-sample-per-symbol sequence directly
        RX_resample_data=zeros(max(ceil(Ndata*P./Q/2)),N_chnn);
        for Chnn=1:N_chnn
            temp= resample(RX_data(:,Chnn), P, 2*Q(Chnn));
            RX_resample_data(1:length(temp),Chnn)= temp;
        end
    end
end


