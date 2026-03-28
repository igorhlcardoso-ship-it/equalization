A = 1; Delta = 1;
M = 4; % QPSK => M = 4
F = 8; Omega = 3; 
Beta = 0.3;
Kc = 4; % number of paths ò RAYLEIGH CHANNEL
Heq = [0.225 0.54 0.675 -0.45]; % impulse response of EDTC
%Heq = [2 -0.4j , 1.5+1.8j, 1, 1.2 -1.3j, 0.8+1.6 j];
fig1 = figure ();
stem(abs(Heq), 'filled');
xlim ([0 10]) ;
title ('Impulse response of equivalent discrete - time channel','FontSize', 12);
xlabel ('Normalized time \times T_{s}','FontSize', 10);
ylabel ('|V| (dB)', 'FontSize', 10); grid on;

fig2 = figure ();
freqz(Heq);
title ('Frequency response of equivalent discrete - time channe', 'FontSize', 12);
xlabel ('Normalized frequency','FontSize', 10);