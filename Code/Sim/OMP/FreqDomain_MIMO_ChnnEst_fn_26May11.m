%-----------------------------------------------------
%  Frequency-Domain Channel Estimation for OFDM
% 
%  Created by Jun Tao
%  require only one OFDM symbol for MIMO channel 
%  estimation
%-----------------------------------------------------

function h = FreqDomain_MIMO_ChnnEst_fn_26May11(y, S, Nt, Nr, Ns, K, L, sc_idx, SNR)
%K: subcarrier number
%Ns: oversampling factor

N= Nt; %N: Tx antenna number
M= Nr; %M: Rx antenna number; 
Y= zeros(M,K);
Np= length(sc_idx); %# of pilot subcarriers

%perform FFT of receive TD OFDM symbols
for m=1:Nr
    %fft in Matlab has no normalization, so scale by 1/sqrt(K)
    Y(m,:)= 1/sqrt(K)*fft(y(m,:)); 
    %scatter(real(Y(m,:)),imag(Y(m,:)),'.')
end

%form training matrix
Fp= zeros(Np,L);
for p= 1: Np
    for l= 1: L
        Fp(p,l)= exp(-j*2*pi*(sc_idx(p)-1)*(l-1)/K)/sqrt(K);
    end
end
Sp= [];
for n= 1: N
    Sp= [Sp diag(S(n,sc_idx))*Fp];
end
Sp= Sp*sqrt(K);

%estimated time-domain MIMO channel coefficients
h= [];
for m= 1: M
    h_temp= pinv(Sp)*Y(m,sc_idx).';
    h= [h; h_temp.'];
end
h= sqrt(Ns)*h;

% h(:,1:L)= h(:,L:-1:1);
% h(:,L+1:2*L)= h(:,2*L:-1:L+1);

return
