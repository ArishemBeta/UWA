%--------------------------------------------------------------------
% MAP convolutional decoding
% Created by Jun Tao at Missouri University of Science and Technology
% Feb 2, 2009
%--------------------------------------------------------------------
% 
% Convolutional Encoder (4,[17,13]) 
% Modified on Feb.27, 2009
%
%--------------------------------------------------------------------

function [LLR_info LLe_cod LL_cod] = MAPConvDecoder_plot(LLa_cod,K);

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

A= [-100*ones(Nstate,1) zeros(Nstate,K)]; %initial value for "forward" vector: log(alpha)
A(1,1)= 0;
B= [zeros(Nstate,K) -100*ones(Nstate,1)]; %initial value for "backward" vector: log(beta)
B(1,K+1)= 0;
C= zeros(Nstate,Nstate,K);

%a priori LL for information bits 
logPa= log(0.5)*ones(2,K);
%LL for coded bits
logPb= LLa_cod;

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
    if(0)
        LLR_cod1(k) = wm1+log(temp3)-zm1-log(temp4);
        LLR_cod2(k) = wm2+log(temp5)-zm2-log(temp6); 
    else
        LL_cod1(1,k) = wm1+log(temp3); %bit 1 equal to 0
        LL_cod1(2,k) = zm1+log(temp4); %bit 1 equal to 1
        c= logsum([LL_cod1(1,k) LL_cod1(2,k)]); LL_cod1(1,k)= LL_cod1(1,k)- c; LL_cod1(2,k)= LL_cod1(2,k)- c;
        LL_cod2(1,k) = wm2+log(temp5); %bit 2 equal to 0
        LL_cod2(2,k) = zm2+log(temp6); %bit 2 equal to 1
        c= logsum([LL_cod2(1,k) LL_cod2(2,k)]); LL_cod2(1,k)= LL_cod2(1,k)- c; LL_cod2(2,k)= LL_cod2(2,k)- c;
    end
end

LL_cod= zeros(2,2*K);
LL_cod(:,1:2:end)= LL_cod1;
LL_cod(:,2:2:end)= LL_cod2;

LLe_cod= LL_cod- LLa_cod;

return

%normalize again
norm1= zeros(1,2*K);
for n= 1: 2*K
    norm1(n)= logsum(LLe_cod(:,n));
end
LLe_cod= LLe_cod-ones(2,1)*norm1;


return

