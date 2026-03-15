p_values = 0:0.01:0.5;
NtxBits = 10000; %we are generating 10000 bits information
N = 3 ; % no. of times we are repeating each bit of information
BER_majority = zeros(size(p_values));
BER_minority = zeros(size(p_values));
         for i = 1:length(p_values)
             p = p_values(i); 
        
             % (b) Generate transmitted bits
             s_tx = randi([0 1], 1, NtxBits);
        
            % (c) Repeat each bit N times
             NcodedBits = N * NtxBits;% no. of bits in the channel
             tx_repeated = repmat(s_tx, N, 1); % make the array tx repeat N * 1 times
             s_tx_coded = reshape(tx_repeated, 1, NcodedBits); % take the first bits of tx_repeated and make another of it
        
            % (d) Generate BSC noise
             bsc_noise = rand(1, NcodedBits) > (1 - p);
        
            % (e) Pass through BSC using XOR
             r = xor(s_tx_coded, bsc_noise);
        
             r_matrix = reshape(r, N, NtxBits);
             ones_count = sum(r_matrix, 1);
             threshold = ceil(N/2);
             s_tx_decoded = sum(r_matrix, 1) < threshold;
        
             % Majority voting (Bayesian / ML)
             s_dec_majority = ones_count >= threshold;
             % Minority voting
             s_dec_minority = ones_count < threshold;

             BER_majority(i) = sum(s_dec_majority ~= s_tx) / NtxBits;
             BER_minority(i) = sum(s_dec_minority ~= s_tx) / NtxBits;
             
        end
     
    plot(p_values, theory_error, 'LineWidth', 2);
    plot(p_values, sim_error, '^', 'MarkerSize', 6);


figure;
plot(p_values, BER_majority, 'b-o', 'LineWidth', 2); hold on;
plot(p_values, BER_minority, 'r-^', 'LineWidth', 2);

xlabel('p');
ylabel('Bit Error Rate');
title('Decoder Comparison for Repeat-by-3 Code over BSC');
legend('Majority Voting (Optimal)', ...
       'Minority Voting');
grid on;