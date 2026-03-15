clc; clear; close all;
K = 4; N = 9;
k = sqrt(K); n = sqrt(N);                               
Nsim = 1e4;
qset  = 0:0.1:1;                                           
psucc = zeros(1, length(qset));

for qi = 1:length(qset)
    q = qset(qi); Tsum = 0;
    for sim = 1:Nsim
        m = double(rand(1,K) < 0.5);
        M       = reshape(m, k, k);
        row_par = mod(sum(M,2), 2);
        M_ext   = [M, row_par];
        col_par = mod(sum(M_ext,1), 2);
        c       = [M_ext; col_par]; c = c(:)';           
        y = double(rand(1,N) < q);
        r = c; r(y==1) = -1;                               
        R = reshape(r, n, n);                            
        for iter = 1:n
            for i = 1:n
                row = R(i,:);
                if sum(row==-1)==1, R(i,row==-1)=mod(sum(row(row~=-1)),2); end
            end
            for j = 1:n
                col = R(:,j);
                if sum(col==-1)==1, R(col==-1,j)=mod(sum(col(col~=-1)),2); end
            end
        end

        m_hat = reshape(R(1:k,1:k),1,K);
        if ~any(m_hat==-1) && isequal(m_hat,m), Tsum=Tsum+1; end
    end
    psucc(qi) = Tsum/Nsim;
end

% For each pattern, check if iterative decoder succeeds, weight by q^e*(1-q)^(N-e)
p_fine       = 0:0.01:1;                                  
psucc_theory = zeros(1, length(p_fine));

% precompute which of the 2^9 = 512 erasure patterns are decodable (independent of q)
decodable = false(1, 2^N);
for pat = 0:2^N-1
    e    = double(dec2bin(pat, N)-'0');                    % erasure pattern, length N
    R_th = reshape(e, n, n);                               % 1=erased, 0=known
    % run iterative decoder on erasure pattern (values don't matter, just positions)
    for iter = 1:n
        for i = 1:n
            row = R_th(i,:);
            if sum(row)==1, R_th(i,row==1)=0; end          % recover erased bit
        end
        for j = 1:n
            col = R_th(:,j);
            if sum(col)==1, R_th(col==1,j)=0; end
        end
    end
    decodable(pat+1) = all(R_th(:)==0);                    % success if no erasures remain
end

% for each q, sum P(pattern) over decodable patterns
for qi = 1:length(p_fine)
    q   = p_fine(qi);
    val = 0;
    for pat = 0:2^N-1
        e_count = sum(double(dec2bin(pat,N)-'0'));          
        val = val + decodable(pat+1) * q^e_count * (1-q)^(N-e_count);
    end
    psucc_theory(qi) = val;
end

figure;
plot(p_fine, psucc_theory, 'g-',  'LineWidth', 2);  hold on; 
plot(qset,   psucc,        'r:',  'LineWidth', 1.5);         
plot(qset,   psucc,        'bs',  'MarkerFaceColor','b','MarkerSize',7);
grid on;
xlabel('p'); ylabel('Prob');
title({'Probability of Successful Decoding for (9,4) Product Code',...
       'Red Dotted Line connecting Blue Squares is Monte Carlo Simulation Result',...
       'Green Solid Line is Theoretical Analysis'});
xlim([0 1]); ylim([0 1.2]);
legend('Theoretical','Simulation','','Location','southwest');