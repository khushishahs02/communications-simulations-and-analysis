p_values = 0:0.01:0.5;
NtxBits = 10000; %we are generating 10000 bits information
N_values = [3 5 7]; % no. of times we are repeating each bit of information
figure; hold on;


for n_idx = 1:length(N_values)

    N = N_values(n_idx);

    sim_error = zeros(size(p_values));
    theory_error = zeros(size(p_values));


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
            threshold = ceil(N/2);
            s_tx_decoded = sum(r_matrix, 1) >= threshold;
        
            num_errors = sum(s_tx ~= s_tx_decoded);
            sim_error(i) = num_errors / NtxBits;
        
            theory_error(i) = 0;
            for k = threshold:N
              theory_error(i) = theory_error(i) + nchoosek(N,k)*p^k*(1-p)^(N-k);
            end
        
         end

    plot(p_values, theory_error, 'LineWidth', 2);
    plot(p_values, sim_error, '^', 'MarkerSize', 6);

end


xlabel('p');
ylabel('Error Probability');
title('Repeat-by-N Coding over BSC');
legend('N=3 Theory','N=3 Sim', ...
       'N=5 Theory','N=5 Sim', ...
       'N=7 Theory','N=7 Sim');
grid on;