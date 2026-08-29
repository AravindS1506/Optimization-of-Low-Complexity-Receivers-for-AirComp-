N=20;
Q=10;
O=2:2:10;
mse_th=[2.17141329133112	0.210333280998076	0.0202744578690764	0.00212906758057531	0.000205353892068602];
mse_sim=[0.390446996772141	0.0229945176275295	0.00132539859081514	9.92962088628929e-05	5.81838142814482e-06];
%mse_th=[3.088885	1.550503553	0.815118311	0.50006458	0.331076579	0.203016714	0.141905432	0.10288828	0.074853651	0.053545757	0.041634144];
%mse_sim=[0.345646297	0.160251947	0.079275497	0.051467766	0.037224043	0.020901918	0.014935959	0.010523017	0.008929156	0.005611152	0.004539107];
%mse_zf=[25.6385019	18.32841293	15.00136586	12.37090975	10.40183229	8.184453923];
%mse_hybrid=[25.56614047	18.92289366	15.78058612	13.16751667	11.40323143	9.198520744];
%mse_zfas=[27.99128198	21.33735867	20.68024317	19.69033614	20.06079388	20.66924435];
figure;
set(gcf, 'Units', 'inches', 'Position', [1, 1, 6, 3.5]);
semilogy(O,mse_th,'-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(O,mse_sim,'-^', 'LineWidth', 2, 'MarkerSize', 8)
%plot(O, mse_zf, '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
%plot(O, mse_hybrid, '-^', 'LineWidth', 2, 'MarkerSize', 8);
%plot(O, mse_zfas, '-v', 'LineWidth', 2, 'MarkerSize', 8);

%xlabel('Number of Antennas (N)');
%ylabel('MSE');
%xlabel('Quantization Levels ($N_q$)', 'Interpreter', 'latex');
xlabel('Number of RF chains (O)');
ylabel('Log Approximation Error');
legend('Theoretical Upper Bound','Simulated Error');
%legend('ZF', 'Hybrid Approximation','ZF+AS');
grid off;

set(gca, 'FontSize', 12);  % axes ticks font size
xticks(O);
grid off;

% Save as EPS
print(gcf, 'mse_vs_rfchains', '-depsc2');