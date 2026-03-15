clc; clear; close all;
q_set_ana = 0:0.01:1;
A_e = [1, 9, 36, 84, 117, 81, 27, 0, 0, 0];

p_succ_ana = zeros(size(q_set_ana));

for i = 1:length(q_set_ana)
    q = q_set_ana(i);
    prob_success = 0;
    for e = 0:9
        prob_success = prob_success + A_e(e+1) * (q^e) * ((1-q)^(9-e));
    end
    p_succ_ana(i) = prob_success;
end

%% Final Plot: Comparing Analytical and Monte Carlo Simulations Results
figure;
hold on; grid on;
% Plot theoretical analysis in a green solid line
plot(q_set_ana, p_succ_ana, 'g-', 'LineWidth', 2, 'DisplayName', 'Theoretical Analysis');
% Assuming p_succ_sim from the N_sim = 1000 run in Section 1.3 is available:
% Plot simulation result with red dotted line connecting blue squares
plot(q_set, p_succ_sim, 'r:s', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', ...
'LineWidth', 1.5, 'DisplayName', 'Monte Carlo Simulation');
title('Probability of Successful Decoding for (9,4) Product Code');
xlabel('Erasure Probability q');
ylabel('Success Probability p_{succ}(q)');
legend('Location', 'northeast');
ylim([0 1.2]);
xlim([0 1]);