a_hat = (randn(4,1) + 1i*randn(4,1)) / sqrt(2);
A_hat = a_hat * a_hat';
[~, min_idx] = mink(diag(A_hat), 3);
disp(sum(diag(A_hat(min_idx,min_idx))))
disp(A_hat)