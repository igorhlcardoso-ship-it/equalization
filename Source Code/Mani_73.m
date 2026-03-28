clc
clear all, close all;

A = 1; Delta = 1;
M = 4; % QPSK => M = 4
K = 2;
F = 8; Omega = 3; 
Beta = 0.3;
Kc = 4; % number of paths ò RAYLEIGH CHANNEL
Heq = [0.225 0.54 0.675 -0.45]; % impulse response of EDTC
%Heq = [2 -0.4j , 1.5+1.8j, 1, 1.2 -1.3j, 0.8+1.6 j];
L = 10000; % number of symbols
Lp = 256; % number of pilot symbols
lambda = 100; % propagation delay , number of null symbols
Es_N0 = 0:1:10; % SNR
H = zeros(1, length(Heq)*F);
H(1:8: end -7) = Heq;
he = filtrerootcos(Beta, F, Omega); % create filter
he_matched = he(end:-1:1);
t0_delay = 2*Omega*F;

%% Generate random PILOT
rng (1237) ;
perm = mod(randperm(K*Lp), 2) > 0;
preambule = Bit2SymbolMappingQPSKGray(A, perm); % pilot symbols
preambule_ech = seqsymbech(preambule, F);
preambule_signal = conv(preambule_ech, he); % pilot signal

Peb = zeros(1, length(Es_N0));

for i=1:length(Es_N0)
    SNR(i) = i-1;
    fprintf(1, 'Es/N0 = %.4f dB\n', SNR(i));
    v = (((A*Delta)^2)/4)*10^(-SNR(i)/10);					%noise variance 
    Number_ErrosBits = 0;
    
    %% Transmiter
    m = rand(1, L*K) < 0.5;
    symbol = Bit2SymbolMappingQPSKGray(A,m); 
    frame = [preambule, symbol]; 
    d = [zeros(1,lambda), frame];
    
    %% Modulation
    d_ech = seqsymbech(d, F);
    trans_signal = conv(d_ech, he)*((Delta*Delta)/2);
    
    %% CHANNEL
    temp = filter(H, 1, trans_signal);
    rev_signal = AWGN(Delta, v, temp);

    %% Matched filtering
    rev_signal = conv(rev_signal, he_matched);

    %% Correlation
    corr = conv(rev_signal,flip(conj(preambule_signal)));
    [peak_max,index_peak] = max(abs(corr));
    tau_est = (index_peak-Lp*F-t0_delay/2+1);

    %% SAMPLING
    z_est = zeros(1, Lp+L);       
    for n_sym = 1:length(z_est)
        z_est(n_sym) = rev_signal(F*(n_sym-1) + tau_est);  %with take into account delay
    end
    
    %% DEMODULATION
    d_est = MLSymbolDetectorQPSKlowCPLX(A, Delta, z_est(Lp+1:end));
    m_est = Symbol2BitsDemappingQPSKGray(A, Delta, d_est);
    
    %check if the decided symbol is right or wrong
    Number_ErrosBits = Number_ErrosBits + sum(xor(m, m_est));
    Peb(i) = Number_ErrosBits/(L*K);

end 

%% Plot BER
figure;
Pb_PSK=0.5*erfc(sqrt(10.^((SNR-3)/10)));%Es/N0 = 2Eb/N0 (QPSK)
semilogy(SNR, Pb_PSK, 'r-', SNR, Peb, 'b--^');
title("Bit Error Probability Curve for QPSK", 'FontSize', 12);
xlabel('Es/No (dB)'); ylabel('BER'); grid on;
legend('Theorique','Selective Frequency Channel');
ylim([1 0.00001]);
