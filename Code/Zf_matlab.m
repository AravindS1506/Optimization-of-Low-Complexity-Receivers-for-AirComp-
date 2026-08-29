clear;
%% Simulation Parameters
eth = 0.0001;
Pmax = 10;
sigma_sqr = 1;
rho = 1.5;
numTrials = 5;       
numNodes = 50;          
sinkPos = [100, 0];
areaSize = 200;             
freq = 2.4e9;
c = 3e8;
L = 1000;
antennaHeight = 1.5;
ref_dist = 90;           % Reference distance where path loss is 80 dB
ref_PL_dB = 80;          % Path loss at reference distance (equivalent to -80 dB gain)
rx_gain_dB = 70;         % Receiver power amplification
N=30;
%% Simulation
N_values = 50:50:300;
mean_MSE = zeros(1, length(N_values));
for idx = 1:length(N_values)
    %N = N_values(idx);
    numNodes=N_values(idx);
    MSE_values = zeros(1, numTrials);
    for i = 1:numTrials
        % Generate channel coefficients
        h = generateChannelCoefficients(N, numNodes, areaSize, sinkPos, ref_dist, ref_PL_dB, rx_gain_dB);
        
        % Run optimization algorithm
        MSE_values(i) = ZF(h, N, numNodes, sigma_sqr, rho, eth, Pmax, L);
    end
    mean_MSE(idx) = mean(MSE_values);
end
figure;
plot(N_values, mean_MSE, '-o', 'LineWidth', 2, 'MarkerSize', 8);
%ylim([6, 17]); % Set y-axis limits

xlabel('Number of Antennas (N)');
ylabel('Mean MSE [log scale]');
title('Mean MSE vs. Number of Antennas');
function A_ret = ZF(hk, N, K, sigma_sqr, rho, eth, Pmax, L) 
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
    A_ret = min(MSE);
end


function MSE = compute_mse(a, b, hk, K, sigma_sqr)
    MSE = 0;
    for k = 1:K
        MSE = MSE + abs(a' * hk(:,k) * b(1,k)-1)^2;
    end
    MSE = MSE + sigma_sqr * norm(a)^2;
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