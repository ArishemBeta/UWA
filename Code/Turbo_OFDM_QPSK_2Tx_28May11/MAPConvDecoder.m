%--------------------------------------------------------------------
% MAP convolutional decoding
% Created by Jun Tao
%--------------------------------------------------------------------
% 
% Convolutional Encoder (4,[17,13]) 
%
%--------------------------------------------------------------------

function [LLR_info LLRe_cod] = MAPConvDecoder(LLRa_cod,K);

%avoiding log of zero
if(0)
    InfMin= 1e-320;
else
    InfMin= 0;    
end

T214 = poly2trellis(4,[17 13]);
tblen = 12;
Nstate= T214.numStates; %# of states of convolutional encoder
nextStates= T214.nextStates; %next state matrix for info. bit
nextStates1= [0 4; 4 0; 5 1; 1 5; 6 2; 2 6; 3 7; 7 3]; %next state matrix for coded bit 1
nextStates2= [0 4; 4 0; 5 1; 1 5; 2 6; 6 2; 7 3; 3 7]; %next state matrix for coded bit 2
preStates= [0 1; 2 3; 4 5; 6 7; 0 1; 2 3; 4 5; 6 7]; %previous state matrix
numInputSymbols= T214.numInputSymbols;
outputs= T214.outputs; %output matrix

A= [-1000*ones(Nstate,1) zeros(Nstate,K)]; %initial value for "forward" vector: log(alpha)
A(1,1)= 0;
B= [zeros(Nstate,K) -1000*ones(Nstate,1)]; %initial value for "backward" vector: log(beta)
B(1,K+1)= 0;
C= zeros(Nstate,Nstate,K);

LLRa_info= zeros(1,K); %initial LLR for information bits 
Pa(1,:)= 1./(1+exp(-LLRa_info));  %prob. of info. bit equaling to 0
%Pa(1,K-tblen+1:K)= 1;  %prob. of info. bit equaling to 0
Pa(2,:)= 1-Pa(1,:);
%Pa(2,K-tblen+1:K)= 1;  %prob. of info. bit equaling to 0
logPa= log(InfMin+Pa);

Pb(1,:)= 1./(1+exp(-LLRa_cod)); %prob. of coded bit equaling to 0
Pb(2,:)= 1-Pb(1,:);
logPb= log(InfMin+Pb);
if(0)
    for k=1: length(LLRa_cod)
        if(Pb(1,k)<exp(-10)) 
            logPb(1,k)= -1000;
            logPb(2,k)= 0;
        elseif(Pb(1,k)>1-exp(-10))
            logPb(1,k)= 0;
            logPb(2,k)= -1000;
        else
            logPb(:,k)= log(Pb(:,k));
        end
    end 
end

%calculate all transition metric at time 1 to K
for k= 1: K
    for s= 1: Nstate
        s0= nextStates(s,1)+1; %input 0
        s1= nextStates(s,2)+1; %input 1
        C(s,s0,k)= logPa(1,k)+logPb(bitget(outputs(s,1),2)+1,2*k-1)+...
                            logPb(bitget(outputs(s,1),1)+1,2*k);
        C(s,s1,k)= logPa(2,k)+logPb(bitget(outputs(s,2),2)+1,2*k-1)+...
                            logPb(bitget(outputs(s,2),1)+1,2*k);
    end
end
    
 %calculate the forward metric at time 2 to K
 for k= 2: K 
    for s= 1: Nstate
        s0= preStates(s,1)+1; 
        s1= preStates(s,2)+1;
        A(s,k)= max(A(s0,k-1)+C(s0,s,k-1),A(s1,k-1)+C(s1,s,k-1))+... 
               log(1+exp(-abs(A(s0,k-1)+C(s0,s,k-1)-A(s1,k-1)-C(s1,s,k-1))));
    end
 end

%calculate the forward metric at time K to 1
for k= K: -1 : 1 
    for s= 1: Nstate
        s0= nextStates(s,1)+1; %input 0
        s1= nextStates(s,2)+1; %input 1
        B(s,k)= max(B(s0,k+1)+C(s,s0,k),B(s1,k+1)+C(s,s1,k))+... 
               log(1+exp(-abs(B(s0,k+1)+C(s,s0,k)-B(s1,k+1)-C(s,s1,k))));
    end
end   

%-------calculate LLR for information bit and coded bit-------
for k =1:K
    for s= 1: Nstate
        w(s)= A(s,k)+C(s,nextStates(s,1)+1,k)+B(nextStates(s,1)+1,k+1);
        w1(s)= A(s,k)+C(s,nextStates1(s,1)+1,k)+B(nextStates1(s,1)+1,k+1);
        w2(s)= A(s,k)+C(s,nextStates2(s,1)+1,k)+B(nextStates2(s,1)+1,k+1);
        
        z(s)= A(s,k)+C(s,nextStates(s,2)+1,k)+B(nextStates(s,2)+1,k+1);
        z1(s)= A(s,k)+C(s,nextStates1(s,2)+1,k)+B(nextStates1(s,2)+1,k+1);
        z2(s)= A(s,k)+C(s,nextStates2(s,2)+1,k)+B(nextStates2(s,2)+1,k+1);
    end
    wm= max(w); wm1= max(w1); wm2= max(w2);
    zm= max(z); zm1= max(z1); zm2= max(z2);
    
    temp1= 0; temp2= 0;
    temp3= 0; temp4= 0;
    temp5= 0; temp6= 0;
    for s= 1: Nstate
        temp1= temp1+ exp(w(s)-wm);
        temp2= temp2+ exp(z(s)-zm);
        temp3= temp3+ exp(w1(s)-wm1);
        temp4= temp4+ exp(z1(s)-zm1);
        temp5= temp5+ exp(w2(s)-wm2);
        temp6= temp6+ exp(z2(s)-zm2);
    end
    LLR_info(k) = wm+log(temp1)-zm-log(temp2);
    LLR_cod1(k) = wm1+log(temp3)-zm1-log(temp4);
    LLR_cod2(k) = wm2+log(temp5)-zm2-log(temp6); 
end

LLR_cod= zeros(1,2*K);
LLR_cod(1:2:end)= LLR_cod1;
LLR_cod(2:2:end)= LLR_cod2;

LLR_cod;
LLRa_cod;
size(LLRa_cod)
size(LLR_cod)
LLRe_cod= LLR_cod- LLRa_cod;

return

