N = 10;
O = 2:2:10;
Q = ceil(pi * 10.^((2 + log10(N)) ./ O));
plot(O, Q, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('O'); ylabel('Q');
grid on;
