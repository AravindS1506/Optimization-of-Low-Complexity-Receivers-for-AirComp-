clear;
%% Simulation Parameters
N = 15;
O=5;
eth = 0.001;
eps=0.001;
Pmax = 10;
sigma_sqr = 1;
rho = 2;
numTrials = 10;       
numNodes = 1;          
sinkPos = [100, 0];
areaSize = 200;             
freq = 2.4e9;
c = 3e8;
L = 500;
antennaHeight = 1.5;
ref_dist = 90;           
ref_PL_dB = 80;          
rx_gain_dB = 70;         

%% Simulation
N_values = 5:5:15;
mean_MSE = zeros(1, length(N_values));
for idx = 1:length(N_values)
    MSE_values = zeros(1, numTrials);
    disp(idx)
    for i = 1:numTrials
        % Generate channel coefficients
        h = generateChannelCoefficients(N_values(idx), numNodes, areaSize, sinkPos, ref_dist, ref_PL_dB, rx_gain_dB);
        %if mod(i,20) == 0

        %end
        % Run optimization algorithm
        MSE_values(i) = ZF(h, N_values(idx), numNodes, sigma_sqr, rho, eth,Pmax, L,O);
    end
    mean_MSE(idx)=mean(MSE_values);
end
% Display average MSE
%disp(['Average MSE: ', num2str(mean(MSE_values))]);
figure;
plot(N_values, mean_MSE, '-o', 'LineWidth', 2, 'MarkerSize', 8);
%set(gca, 'YScale', 'log'); % Set y-axis to log scale
%ylim([0, 6]); % Set y-axis limits

xlabel('Number of Antennas (N)');
ylabel('Mean MSE');
title('Mean MSE vs. Number of Antennas');
grid on;
function A_ret = ZF(hk, N, K, sigma_sqr, rho, eth, Pmax, L,O) 
    % Initialize
    a_hat = rand(N,1) + 1i*rand(N,1);
    A_hat = a_hat * a_hat';
    
    % Precompute Hk matrices
    for k = 1:K
        Hk(:,:,k) = hk(:,k) * hk(:,k)';
    end
    
    % Iterative optimization using DC approach
    MSE_prev = inf;
    while true
        % Compute subgradient
        [U, ~, ~] = svd(A_hat);
        v = U(:,1);
        subgrad = v * v';

        % Solve the convex subproblem for ZF (all signals aligned)
        cvx_begin quiet
            variable A(N,N) semidefinite
            minimize(sigma_sqr*trace(A) + rho*trace(A*(eye(N)-subgrad)))
            for k = 1:K
                    real(trace(A * Hk(:,:,k))) >= 1 / Pmax;
            end
        cvx_end 
        
        % Check convergence
        MSE_curr = cvx_optval;
        if abs(MSE_curr - MSE_prev)/abs(MSE_curr) <= eth 
            break
        end
        
        MSE_prev = MSE_curr;
        %[U,S,~]=svd(A);
        %U1 = U(:,1);
        %sigma1 = S(1,1);
        %ahat_t = U1 * sqrt(sigma1);
        %A_hat=ahat_t*ahat_t';
        A_hat=A;
    end

    % Gaussian randomization
    [U, Sigma] = svd(A);
    a = zeros(N, L);
    MSE = zeros(1, L);
    
    for j = 1:L
        % Generate random vector
        n = (randn(N,1) + 1i*randn(N,1))/sqrt(2);
        a(:,j) = U * sqrt(Sigma) * n;
        
        % Normalize to satisfy ZF constraints
        min_val = min(abs(a(:,j)' * hk));
        a(:,j) = (a(:,j) / min_val) * (1 / sqrt(Pmax));
        
        % Compute b_k for each node to achieve alignment
        for k = 1:K
            b(j,k) = 1/(a(:,j)' * hk(:,k));
        end
        
        % For ZF, MSE is just the noise term since all signals are aligned
        MSE(j) = sigma_sqr * norm(a(:,j))^2;
    end
    
    % Find best result
    [A_ret,idx] = min(MSE);
    [~,~,a_hybrid]=unit_modulus_factorization(a(:,idx),O);
    min_val = min(abs(a_hybrid' * hk));
    a_hybrid = (a_hybrid / min_val) * (1 / sqrt(Pmax));
    %disp(norm(a_hybrid-a(:,idx))^2)
    A_ret=sigma_sqr*norm(a_hybrid)^2;
end

%{
function [W, ao, an] = unit_modulus_factorization(v, O)
% UNIT_MODULUS_FACTORIZATION Factorizes vector v ≈ W*ao with |W(i,j)| = 1
% Uses coordinate descent for improved convergence

N = length(v);
max_iter = 1000;  % Maximum number of iterations
tol = 1e-6;       % Tolerance for convergence
lambda = 1e-8;    % Regularization parameter

% Multiple random restarts
num_restarts = 3;
best_residual = inf;
W_best = exp(1i * 2*pi * rand(N, O));
ao_best = zeros(O,1);

for restart = 1:num_restarts
    W = exp(1i * 2*pi * rand(N, O)); 
    residuals = zeros(max_iter, 1);
    
    for iter = 1:max_iter
        % Compute ao using regularized least squares
        ao = (W'*W + lambda*eye(O)) \ (W'*v);
        
        % Coordinate descent update for W
        for j = 1:N
            s_j = W(j,:) * ao;  % Current approximation for v(j)
            for k = 1:O
                W_old = W(j,k);
                r_j = v(j) - (s_j - W_old * ao(k));  % Residual excluding W(j,k)*ao(k)
                if abs(ao(k)) > 1e-10  % Avoid issues with zero coefficients
                    W_new = exp(1i * (angle(r_j) - angle(ao(k))));
                else
                    W_new = W_old;  % If ao(k) = 0, keep W(j,k) unchanged
                end
                s_j = s_j - W_old * ao(k) + W_new * ao(k);  % Update s_j incrementally
                W(j,k) = W_new;
            end
        end
        
        % Compute residual
        residual = norm(v - W*ao);
        residuals(iter) = residual;
        
        % Early stopping if improvement is negligible
        if iter > 50 && abs(residuals(iter-50) - residual) < 1e-3*tol
            break
        end
    end
    
    % Track the best solution across restarts
    if residual < best_residual
        best_residual = residual;
        W_best = W;
        ao_best = ao;
    end
end

% Return the best solution found
W = W_best;
ao = (W'*W + lambda*eye(O)) \ (W'*v);  % Recompute ao for consistency
an = W*ao;

end
%}
function [W, ao, an] = unit_modulus_factorization(v, O)
% UNIT_MODULUS_FACTORIZATION Factorizes vector v ≈ W*ao with W(i,j) in Q-th roots of unity
% W entries are constrained to Q evenly spaced unit modulus values

N = length(v);
max_iter = 1000;  % Maximum number of iterations
tol = 1e-6;       % Tolerance for convergence
lambda = 1e-8;    % Regularization parameter
Q=10;
% Define Q possible unit modulus values
theta_Q = 2*pi*(0:(Q-1))/Q;  % Angles: 0, 2*pi/Q, ..., 2*pi*(Q-1)/Q
W_values = exp(1i * theta_Q); % Possible values for W entries

% Multiple random restarts
num_restarts = 1;
best_residual = inf;
W_best = W_values(randi(Q, N, O));  % Initialize W with random Q-th roots
ao_best = zeros(O,1);

for restart = 1:num_restarts
    W = W_values(randi(Q, N, O));  % Random initialization from W_values
    W = reshape(W, N, O);
    residuals = zeros(max_iter, 1);
    
    for iter = 1:max_iter
        % Compute ao using regularized least squares
        ao = (W'*W + lambda*eye(O)) \ (W'*v);
        
        % Coordinate descent update for W
        for j = 1:N
            s_j = W(j,:) * ao;  % Current approximation for v(j)
            for k = 1:O
                W_old = W(j,k);
                r_j = v(j) - (s_j - W_old * ao(k));  % Residual excluding W(j,k)*ao(k)
                if abs(ao(k)) > 1e-10  % Avoid issues with zero coefficients
                    % Compute ideal phase
                    ideal_phase = angle(r_j) - angle(ao(k));
                    % Snap to nearest allowed phase
                    [~, idx] = min(abs(mod(ideal_phase - theta_Q + pi, 2*pi) - pi));
                    W_new = W_values(idx);
                else
                    W_new = W_old;  % If ao(k) = 0, keep W(j,k) unchanged
                end
                s_j = s_j - W_old * ao(k) + W_new * ao(k);  % Update s_j incrementally
                W(j,k) = W_new;
            end
        end
        
        % Compute residual
        residual = norm(v - W*ao);
        residuals(iter) = residual;
        
        % Early stopping if improvement is negligible
        if iter > 50 && abs(residuals(iter-50) - residual) < 1e-3*tol
            break
        end
    end
    
    % Track the best solution across restarts
    if residual < best_residual
        best_residual = residual;
        W_best = W;
        ao_best = ao;
    end
    %disp(W_best)
end

% Return the best solution found
W = W_best;
ao = (W'*W + lambda*eye(O)) \ (W'*v);  % Recompute ao for consistency
an = W*ao;
end


function h = generateChannelCoefficients(N, numNodes, areaSize, sinkPos, ref_dist, ref_PL_dB, rx_gain_dB)
    
    % 1. Node Distribution
    nodes = rand(numNodes, 2) * areaSize;  % Uniform random deployment
    
    % 2. Distance Calculations
    dx = nodes(:,1) - sinkPos(1);
    dy = nodes(:,2) - sinkPos(2);
    distances = sqrt(dx.^2 + dy.^2);
    
    
    % Free-space model (d <= 90m)
    PL_dB = 20*log10(distances/ref_dist) + ref_PL_dB;
    
    
    % 5. Receiver Gain Application
    effective_PL_dB = PL_dB - rx_gain_dB;
    
    % 6. Channel Gain Conversion
    channel_gain = 10.^(-effective_PL_dB/10);  % Linear scale
    
    % 7. Rayleigh Fading Implementation
    h = (randn(N, numNodes) + 1i*randn(N, numNodes))/sqrt(2);  % CN(0,1)
    h = h .* sqrt(reshape(channel_gain, [1, numNodes]));  % Apply path loss scaling
end
