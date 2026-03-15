p = 0:0.01:1;

H_X_given_Y = zeros(size(p));
I_XY = zeros(size(p));
P_error = zeros(size(p));
lambda_all = zeros(4, length(p));
q_all = zeros(4, length(p));


for i = 1:length(p)

    p_now = p(i);

    if p_now == 0
        L = Inf;
    elseif p_now == 1
        L = 0;
    else
        L = (1 - p_now) / p_now;
    end

    H_total = 0; 

    for k = 0:3

        lambda = L^(2*k - 3);
        lambda_all(k+1, i) = lambda;


        if isinf(lambda)
            q_tilde = 1;
        elseif lambda == 0
            q_tilde = 0;
        else
            q_tilde = lambda / (1 + lambda);
        end

        q_all(k+1, i) = q_tilde;
        if (q_tilde == 0 || q_tilde == 1)
            H_cond = 0;
        else
            H_cond = -q_tilde*log2(q_tilde) - (1 - q_tilde)*log2(1 - q_tilde);
        end

        P_Y = 0.5*((1 - p_now)^k * p_now^(3 - k) + p_now^k * (1 - p_now)^(3 - k));

        H_total = H_total + nchoosek(3, k) * P_Y * H_cond;

    end 

    H_X_given_Y(i) = H_total;

    I_XY(i) = 1 - H_X_given_Y(i);

    P_error(i) = binopdf(2, 3, p_now) + ...
                 binopdf(3, 3, p_now);

end

figure;
plot(p, I_XY, 'LineWidth',2);
xlabel('p');
ylabel('I(X;Y)');
title('Mutual Information (N=3)');
grid on;

figure;
plot(p, H_X_given_Y, 'LineWidth',2);
xlabel('p');
ylabel('H(X|Y)');
title('Conditional Entropy (N=3)');
grid on;

figure;
plot(p, P_error, 'LineWidth',2);
xlabel('p');
ylabel('Decoder Error Probability');
title('Decoder Failure Probability (N=3)');
grid on;

figure;
semilogy(p, lambda_all(1,:), 'LineWidth',2); hold on;
semilogy(p, lambda_all(2,:), 'LineWidth',2);
semilogy(p, lambda_all(3,:), 'LineWidth',2);
semilogy(p, lambda_all(4,:), 'LineWidth',2);

xlabel('p');
ylabel('\lambda^{(N=3)}');
title('Likelihood Ratio for N=3');
legend('k=0','k=1','k=2','k=3');
grid on;

figure;
plot(p, q_all(1,:), 'LineWidth',2); hold on;
plot(p, q_all(2,:), 'LineWidth',2);
plot(p, q_all(3,:), 'LineWidth',2);
plot(p, q_all(4,:), 'LineWidth',2);

xlabel('p');
ylabel('q_tilde^{(N=3)}');
title('Posterior Probability for N=3');
legend('k=0','k=1','k=2','k=3');
grid on;
