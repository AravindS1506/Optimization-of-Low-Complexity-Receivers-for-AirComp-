clear;
%% Simulation Parameters
N = 15;
O=5;
eth = 0.001;
eps=0.001;
Pmax = 10;
sigma_sqr = 1;
rho = 2;
numTrials = 1;       
numNodes = 50;          
sinkPos = [100, 0];
areaSize = 200;             
freq = 2.4e9;
c = 3e8;
L = 500;
antennaHeight = 1.5;
ref_dist = 90;           
ref_PL_dB = 80;          
rx_gain_dB = 70;         
eps=10^-2;
N_values = 5:5:30;
mean_MSE = zeros(1, length(N_values));
mean_MSE_th = zeros(1, length(N_values));
mean_MSE_opt = zeros(1, length(N_values));
for idx = 1:length(N_values)
    MSE_values = zeros(1, numTrials);
    MSE_theory = zeros(1, numTrials);
    MSE_opt = zeros(1, numTrials);
    disp(idx)
    Q=ceil(pi*10^((-log10(eps)+log10(N))/O));
    for i = 1:numTrials
        % Generate channel coefficients
        h = generateChannelCoefficients(N_values(idx), numNodes, areaSize, sinkPos, ref_dist, ref_PL_dB, rx_gain_dB);
        [MSE_values(i),MSE_opt(i),MSE_theory(i)] = ZF(h, N_values(idx), numNodes, sigma_sqr, rho, eth,Pmax, L,O,Q);
    end
    mean_MSE(idx)=mean(MSE_values);
    mean_MSE_th(idx)=mean(MSE_theory);
    mean_MSE_opt(idx)=mean(MSE_opt);
end
% Display average MSE
%disp(['Average MSE: ', num2str(mean(MSE_values))]);
figure;
plot(N_values, mean_MSE, '-o', 'LineWidth', 2, 'MarkerSize', 8);hold on;
%plot(N_values, mean_MSE_th, '-^', 'LineWidth', 2, 'MarkerSize', 8);
plot(N_values, mean_MSE_opt, '-v', 'LineWidth', 2, 'MarkerSize', 8);hold off;
%set(gca, 'YScale', 'log'); % Set y-axis to log scale
%ylim([0, 6]); % Set y-axis limits

xlabel('Number of Antennas (N)');
ylabel('Mean MSE');
title('Mean MSE vs. Number of Antennas');
grid on;

function [MSE_ZF,MSE_opt,MSE_th] = ZF(hk, N, K, sigma_sqr, rho, eth, Pmax, L,O,Q) 
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
        
        
        % For ZF, MSE is just the noise term since all signals are aligned
        MSE(j) = sigma_sqr * norm(a(:,j))^2;
    end
    
    % Find best result
    [MSE_opt,idx] = min(MSE);
    a_min=a(:,idx);
    Q_mat=-sigma_sqr*eye(N);
    for k = 1:K
        inner_prod = a_min' * hk(:,k); % Compute a^H h_k
        b(k) = 1 / inner_prod;
        Q_mat=Q_mat+abs(b(k))^2*(hk(:,k)*hk(:,k)');
    end
    [W,~,a_hybrid]=unit_modulus_factorization(a_min,O,Q,Q_mat);
    MSE_ZF=compute_mse(a_hybrid,b,hk,K,sigma_sqr);
    %lambda_max = eigs(Q_mat, 1, 'largestreal');
    vals = eigs(Q_mat,N, 'largestreal');  % returns top O+1 eigenvalues (unsorted)
    vals = sort(vals, 'descend');
    if N~=O
        %disp(rank(W))
        lambda_max = vals(1);
        MSE_th = norm(a_hybrid - a_min)^2 * (real(lambda_max)) + sigma_sqr*norm(a_min)^2;
    else
        MSE_th=sigma_sqr*norm(a_min)^2;
    end
end


function MSE = compute_mse(a, b, hk, K, sigma_sqr)
    MSE = 0;
    for k = 1:K
        MSE = MSE + abs(a' * hk(:,k) * b(k)-1)^2;
    end
    MSE = MSE+sigma_sqr * norm(a)^2;
end

function [W, ao, an] = unit_modulus_factorization(v, O, Q, H)
% UNIT_MODULUS_FACTORIZATION with weighted norm using matrix H

N = length(v);
max_iter = 1000;
tol = 1e-6;
lambda = 1e-8;

% Define Q possible unit modulus values
theta_Q = 2*pi*(0:(Q-1))/Q;
W_values = exp(1i * theta_Q);

num_restarts = 3;
best_residual = inf;
W_best = W_values(randi(Q, N, O));
ao_best = zeros(O, 1);

for restart = 1:num_restarts
    W = W_values(randi(Q, N, O));
    W = reshape(W, N, O);
    residuals = zeros(max_iter, 1);
    
    for iter = 1:max_iter
        % Weighted least squares using H
        Wh = W' * H;
        ao = (Wh * W + lambda * eye(O)) \ (Wh * v);

        % Coordinate descent update for W
        for j = 1:N
            s_j = W(j,:) * ao;
            for k = 1:O
                W_old = W(j,k);
                r_j = v(j) - (s_j - W_old * ao(k));
                if abs(ao(k)) > 1e-10
                    ideal_phase = angle(r_j) - angle(ao(k));
                    [~, idx] = min(abs(mod(ideal_phase - theta_Q + pi, 2*pi) - pi));
                    W_new = W_values(idx);
                else
                    W_new = W_old;
                end
                s_j = s_j - W_old * ao(k) + W_new * ao(k);
                W(j,k) = W_new;
            end
        end

        % Weighted residual computation
        eps = v - W * ao;
        residual = real(eps' * H * eps);
        residuals(iter) = residual;

        if iter > 50 && abs(residuals(iter - 50) - residual) < 1e-3 * tol
            break
        end
    end

    if residual < best_residual
        best_residual = residual;
        W_best = W;
        ao_best = ao;
    end
end

W = W_best;
Wh = W' * H;
ao = (Wh * W + lambda * eye(O)) \ (Wh * v);
an = W * ao;

end

%{
function [W, ao, an] = unit_modulus_factorization(v, O,Q)
N = length(v);
theta_Q = 2*pi*(0:Q-1)/Q;      % Allowed angles
W_values = exp(1i * theta_Q);  % Allowed complex values on unit circle

W = zeros(N, O);    % Store selected vectors
ao = zeros(O, 1);   % Coefficients
an=zeros(N,1);
residual = v;       % Initial residual

for o = 1:O
    w_best = zeros(N, 1);

    % Greedy selection: align each entry to the closest quantized phase of residual
    for i = 1:N
        phase_i = angle(residual(i));
        [~, idx] = min(abs(wrapToPi(theta_Q - phase_i)));  % Nearest quantized phase
        w_best(i) = W_values(idx);
    end

    % Compute optimal ao for this w
    ao(o) = (w_best' * residual) / N;

    % Update W and residual
    W(:, o) = w_best;
    residual = residual - ao(o) * w_best;
    an=an+ao(o)*w_best;
end
end

%}
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
