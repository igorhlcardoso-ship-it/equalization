A = 1; Delta = 1;
M = 4; % QPSK => M = 4
K = 2;
F = 8; Omega = 3; 
Beta = 0.3;
Kc = 4; % number of paths ò RAYLEIGH CHANNEL
Heq = [0.225 0.54 0.675 -0.45]; % impulse response of EDTC
%Heq = [2 -0.4j , 1.5+1.8j, 1, 1.2 -1.3j, 0.8+1.6 j];
L = 5000; % number of symbols
Lp = 5000; % number of pilot symbols
lambda = 100; % propagation delay , number of null symbols
Es_N0 = 15; % SNR
H = zeros(1, length(Heq)*F);
H(1:8: end -7) = Heq;
he = filtrerootcos(Beta, F, Omega); % create filter
he_matched = he(end:-1:1);
t0_delay = 2*Omega*F;
mu = 0.02;
delta = 15; 
Ng = 30;
gamma_RLS = 0.99;

%% Generate random PILOT
rng (1237) ;
perm = mod(randperm(K*Lp), 2) > 0;
preambule = Bit2SymbolMappingQPSKGray(A, perm); % pilot symbols
preambule_ech = seqsymbech(preambule, F);
preambule_signal = conv(preambule_ech, he); % pilot signal

Peb = zeros(1, length(Es_N0));

for i=1:length(Es_N0)
    SNR(i) = Es_N0;
    fprintf(1, 'Es/N0 = %.4f dB\n', SNR(i));
    v = (((A*Delta)^2)/4)*sum(abs(Heq).^2*10.^(-SNR(i)/10));					%noise variance 
    Number_ErrosBits = 0;
    
    %% Transmiter
    m = rand(1, L*K) < 0.5;
    symbol = Bit2SymbolMappingQPSKGray(A,m); 
    frame = [preambule, symbol]; 
    d = [zeros(1,lambda), frame];
    
    %% Modulation
    d_ech = seqsymbech(d, F);
    trans_signal = conv(d_ech, he);
    
    %% CHANNEL
    temp = filter(H, 1, trans_signal);
    rev_signal = AWGN(Delta, v, temp);

    %% Matched filtering
    rev_signal = conv(rev_signal, he_matched);

    %% Correlation
    corr = conv(rev_signal,flip(conj(preambule_signal)));
    [peak_max,index_peak] = max(abs(corr));
    tau_est = (index_peak-Lp*F-t0_delay/2 +1);

    %% SAMPLING
    z_est = zeros(1, Lp+L);       
    for n_sym = 1:length(z_est)
        z_est(n_sym) = rev_signal(F*(n_sym-1) + tau_est);  %with take into account delay
    end

    %% EQUALIZATION 
    y = zeros(1, Lp);
    zl = zeros(1, Ng);
    g = zeros(1, Ng);
    MSE = zeros(1, Ng);

    for k = Ng:Lp+delta
         zl(end:-1:1) = z_est(k-Ng+1:k);
         y(k) = zl*transpose(g);
         g = g - mu*(y(k) - preambule(k-delta))*conj(zl);
         MSE(k) = gamma_RLS*MSE(k-1) + (1-gamma_RLS)*(abs(y(k) - preambule(k-delta)))^2;
    end

    fig1 = figure();
    scatter(real(z_est), imag(z_est), 10, 'filled'); grid on;
    title('Constellation before equalization', 'FontSize', 12);
    xlabel('I');
    ylabel('Q');
    axis equal;

    fig2 = figure();
    scatter(real(y(end-1000:end)), imag(y(end-1000:end)), 10, 'filled'); grid on;
    title('Constellation after equalization', 'FontSize', 12);
    xlabel('I');
    ylabel('Q');
    axis equal;
    
    fig3 = figure();
    semilogy(MSE);
    grid on; 
    title('Evolution of MSE', 'FontSize', 12);
    legend('\mu = 0.003');
    xlabel('Samples');
 end 