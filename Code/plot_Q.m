N = 10;
O = 2:2:10;  % O = [2, 4, 6, 8, 10]
disp('O:');
disp(O);

exponent = (2 + log10(N)) ./ O;
disp('Exponent:');
disp(exponent);

Q = ceil(pi * 10.^exponent);
disp('Q:');
disp(Q);

plot(O, Q, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('O');
ylabel('Q');
grid on;
