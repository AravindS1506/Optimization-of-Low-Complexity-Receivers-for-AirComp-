clear;
%% Simulation Parameters

num_S=0;
eth = 0.01;
Pmax = 10;
sigma_sqr = 1;
rho = 1;
numTrials = 100;       
numNodes = 50;          
sinkPos = [100, 0];
areaSize = 200;             
freq = 2.4e9;
c = 3e8;
L = 500;
antennaHeight = 1.5;
ref_dist = 90;           % Reference distance where path loss is 80 dB
ref_PL_dB = 80;          % Path loss at reference distance (equivalent to -80 dB gain)
rx_gain_dB = 70;         % Receiver power amplification

%% Simulation
N_values = 5:5:30;
mean_MSE = zeros(1, length(N_values));
for idx = 1:length(N_values)
    disp(idx)
    N = N_values(idx);
    MSE_values = zeros(1, numTrials);
    for i = 1:numTrials
        if(mod(i,5)==0)
            disp(i);
        end
        % Generate channel coefficients
        h = generateChannelCoefficients(N, numNodes, areaSize, sinkPos, ref_dist, ref_PL_dB, rx_gain_dB);
        
        % Run optimization algorithm
        MSE_values(i) = miso(h, N, numNodes, sigma_sqr, rho, eth, Pmax, L);
    end
    mean_MSE(idx) = mean(MSE_values);
end
figure;
plot(N_values, mean_MSE, '-o', 'LineWidth', 2, 'MarkerSize', 8);
%ylim([6, 17]); % Set y-axis limits

xlabel('Number of Antennas (N)');
ylabel('Mean MSE [log scale]');
title('Mean MSE vs. Number of Antennas');
function A_ret = miso(hk, N, K, sigma_sqr, rho, eth, Pmax, L) 
    % Initialize
    num_S=0;
    a_hat_prev = (randn(N,1) + 1i*randn(N,1)) / sqrt(2);
    A_hat = a_hat_prev * a_hat_prev';
    A_t_prev=A_hat;
    a_t_prev=a_hat_prev;
    Hk = zeros(N, N, K);
    % Precompute Hk matrices
    for k = 1:K
        Hk(:,:,k) = hk(:,k) * hk(:,k)';
    end
    A_ret=0;
    % Iterative optimization using DC approach
    while true

        scores=abs(a_t_prev'*hk);
        [~,top_sorted]=sort(scores,'ascend');
        top_sorted=top_sorted(1:num_S);
        not_S = setdiff(1:K, top_sorted);
        % Compute subgradient
        [U, ~, ~] = svd(A_t_prev);
        v = U(:,1);
        subgrad = v * v';

        % Solve the convex subproblem for ZF (all signals aligned)
        cvx_begin quiet
            variable A(N,N) semidefinite hermitian
            miso_tr = 0;
            if ~isempty(top_sorted)
                for i = 1:length(top_sorted)
                    miso_tr = miso_tr + real(trace(A * Hk(:,:,top_sorted(i))));
                end
            end
            minimize(sigma_sqr * trace(A) + rho * trace(A * (eye(N) - subgrad))+ num_S - Pmax * miso_tr)
            subject to
                % Constraints for k in S
                for k = 1:length(top_sorted)
                    real(trace(A * Hk(:,:,top_sorted(k)))) <= 1 / Pmax;
                end
                % Constraints for k not in S
                if ~isempty(not_S)
                    for k = 1:length(not_S)
                        real(trace(A * Hk(:,:,not_S(k)))) >= 1 / Pmax;
                    end
                end
        cvx_end

        %disp(cvx_optval)
        
        % Check convergence
        [U, D] = eig(A');  % Eigen decomposition of A^t

        % Extract the largest eigenvalue and its corresponding eigenvector
        [~, idx] = max(diag(D));
        sigma1 = D(idx, idx);
        U1 = U(:, idx);
        a_hat_t = U1 * sqrt(sigma1);
        %disp(cvx_optval);
        if norm(a_hat_t - a_hat_prev)/norm(a_hat_t) <= eth 
            %disp("Converge");
            break
        end
        a_hat_prev=a_hat_t;
        a = zeros(N, L);
        b = zeros(L, K);
        MSE = zeros(1, L);
        
        for j = 1:L
            % Generate random vector
            n = (randn(N,1) + 1i*randn(N,1))/sqrt(2);
            a(:,j) = U * sqrt(D) * n;
            
            % Normalize to satisfy ZF constraints
            remaining_nodes = setdiff(1:K, top_sorted);
            if ~isempty(remaining_nodes)
                min_val = min(abs(a(:,j)' * hk(:, remaining_nodes)));
                a(:,j) = (a(:,j) / min_val) * (1 / sqrt(Pmax));
            end
            
            % Compute b_k for each node to achieve alignment
            for k = 1:length(top_sorted)
                MSE(j)=MSE(j)+(1-abs(a(:,j)' * hk(:,top_sorted(k))*sqrt(Pmax)))^2;
            end
            
            % For ZF, MSE is just the noise term since all signals are aligned
            MSE(j) = MSE(j)+sigma_sqr * norm(a(:,j))^2;
        end
        
        % Find best result
        [A_ret, min_index] = min(MSE);
        a_t_prev=a(:,min_index);
        A_t_prev = A;        
        num_S=num_S+1;
        if(num_S==K+1)
            break;
        end
    end
end


function MSE = compute_mse(a, b, hk, K, sigma_sqr)
    MSE = 0;
    for k = 1:K
        MSE = MSE + abs(a' * hk(:,k) * b(k)-1)^2;
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