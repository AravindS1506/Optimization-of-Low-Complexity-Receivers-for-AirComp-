data = [
50	100	150	200	250	300
6.035902197	15.08967491	26.1359925	38.26670124	50.68226152	67.98316183
5.946077593	14.88433489	25.7877113	37.67763002	49.91261199	67.0487036
];

K = data(1, :);        % x-axis: node counts
MSE_1 = data(2, :);    % First method's MSE
MSE_2 = data(3, :);    % Second method's MSE

figure;
plot(K, MSE_1, '-o', 'LineWidth', 2); hold on;
plot(K, MSE_2, '-^', 'LineWidth', 2); hold off;
xlabel('Number of Nodes (K)');
ylabel('Mean MSE');
legend('Proposed Algorithm', 'Ideal Case');
grid on;
