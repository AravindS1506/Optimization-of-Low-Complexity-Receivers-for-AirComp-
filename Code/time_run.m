O=2:2:10;
time_zf=[5.968543783	5.849813437	5.878325996	6.098526723	5.834788009];
time_zfas=[10.06909235	8.831538783	7.713535383	7.12003781	5.540798308];

figure;
set(gcf, 'Units', 'inches', 'Position', [1, 1, 6, 3.5]);
plot(O, time_zfas, '-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(O, time_zf, '-s', 'LineWidth', 1.5, 'MarkerSize', 8);

xlabel("Number of RF chains O");
ylabel("Time (s)");

% Set y-axis from 0 to 12, tick every 2
ylim([0 12]);
yticks(0:2:12);

% Set x-axis from 0 to 10 and label according to O values
set(gca, 'FontSize', 12);
xlim([2 10]);
xticks(O);
xticklabels(string(O));

% Add legend
legend("ZF+AS algorithm", "ZF+Our Algorithm", 'Location', 'northeast');

grid on;

