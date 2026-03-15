clearvars; close all;
s = 1;
sigma2 = 1; sigma = sqrt(sigma2);
r_decision_boundary = -1;
Nsim = 10000;
Ncorr = 0;

gamma_dB = 0:1:10; % SNR Range in dB

BER_sim = zeros(size(gamma_dB));
BER_theory = zeros(size(gamma_dB));

for i = 1:length(gamma_dB)
    gamma = 10^(gamma_dB(i)/10);
    
    sigma2 = s^2/(2*gamma);
    sigma = sqrt(sigma2);

    Nerr = 0;
    
    % Monte-Carlo Experimental Trials
    for ksim = 1:Nsim
        Xtx = randi(2,1) - 1; % random generation of X = 0 and X = 1
        if Xtx == 1 % if X = 1, transmit s, else transmit -s volts
            s_tx = s;
        else
            s_tx = -s;
        end
        % Matlab's function randn generates Gaussian distributed random
        % variable with variance of 1. Multiply by sigma to make the variance
        % sigma^2
        n = sigma*randn;
        % received signal
        r = s_tx + n;
        % Bayesian receiver
        if r > r_decision_boundary
            Xhat_rx = 1;
        else
            Xhat_rx = 0;
        end
% Update the counters
        if Xtx ~= Xhat_rx % if the receiver has decided correctly
           Nerr = Nerr + 1; %increment the count of correct bit decision at the receiver
        end
    end
    BER_sim(i) = Nerr / Nsim;
    BER_theory(i) = 0.5 * erfc(sqrt(gamma));

end

% Plot BER vs SNR
semilogy(gamma_dB, BER_sim, 'o-', 'LineWidth', 1.5); hold on;
semilogy(gamma_dB, BER_theory, '--', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('Simulated BER', 'Theoretical BER');
title('BER Performance of BPSK over AWGN Channel');