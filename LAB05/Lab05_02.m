clc; clear; close all;
K = 4; N = 9;                                              
k = sqrt(K); n = sqrt(N);                                  
q = 0.3;                                                  

% Q1: Generate K message bits m ~ Bern(0.5)
m = double(rand(1,K) < 0.5);                             
fprintf('Message bits m = '); disp(m);

% Q2: Encoder - build N=9 encoded bits array c
% Arrange K=4 info bits in k x k (2x2) grid, add parity row and column
M = reshape(m, k, k);                                      % 2x2 info bit matrix
% parity: XOR along rows and columns (even parity)
row_par = mod(sum(M, 2), 2);                               % 2x1 row parities
M_ext   = [M, row_par];                                    % 2x3 extended with row parity
col_par = mod(sum(M_ext, 1), 2);                           % 1x3 col parities
C       = [M_ext; col_par];                                % 3x3 codeword matrix
c       = C(:)';                                           % N=9 encoded bits, row-major
fprintf('Codeword c = '); disp(c);

% Q3: BEC channel - generate N erasure flags y ~ Bern(q)
y = double(rand(1,N) < q);                                 % 1=erased, 0=received
fprintf('Erasure flags y = '); disp(y);

% Q4: Apply erasures to get received word r (-1 denotes erased bit)
r = c;
r(y == 1) = -1;                                            % -1 marks erased positions
fprintf('Received r = '); disp(r);

% Q5: Iterative row-column decoder
R = reshape(r, n, n);                                      % reshape to 3x3 matrix

% iterate up to n times (enough for any recoverable pattern)
for iter = 1:n
    % row-by-row: if exactly one -1 in a row, recover it by XOR of known bits
    for i = 1:n
        row = R(i,:);
        if sum(row == -1) == 1                             % exactly one erasure in row
            R(i, row == -1) = mod(sum(row(row ~= -1)), 2); % recover by even parity
        end
    end
    % column-by-column: same logic
    for j = 1:n
        col = R(:,j);
        if sum(col == -1) == 1                             % exactly one erasure in col
            R(col == -1, j) = mod(sum(col(col ~= -1)), 2);
        end
    end
end

% extract decoded message bits from top-left k x k submatrix
m_hat = reshape(R(1:k, 1:k), 1, K);                       % decoded message bits
fprintf('Decoded m_hat = '); disp(m_hat);

% Q6: Decoder success flag T
if any(m_hat == -1) || ~isequal(m_hat, m)
    T = 0;                                                 % failure: erasures remain or mismatch
else
    T = 1;                                                 % success: decoded matches transmitted
end
fprintf('Decoder success flag T = %d\n', T);