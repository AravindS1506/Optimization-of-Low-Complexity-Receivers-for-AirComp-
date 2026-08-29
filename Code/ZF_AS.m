clear;
%% Simulation Parameters
N = 20;
O=4;
eth = 0.001;
eps=0.001;
Pmax = 10;
sigma_sqr = 1;
rho = 2;
numTrials = 100;       
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

%% Simulation
N_values = 2:2:10;
mean_MSE = zeros(1, length(N_values));
mean_time=zeros(1,length(N_values));
for idx = 1:length(N_values)
    MSE_values = zeros(1, numTrials);
    time=zeros(1,numTrials);
    for i = 1:numTrials
        O=N_values(idx);
        % Generate channel coefficients
        h = generateChannelCoefficients(N, numNodes, areaSize, sinkPos, ref_dist, ref_PL_dB, rx_gain_dB);
        if mod(i, 20) == 0
            fprintf('i = %d, idx = %d\n', i, idx);
        end
        t_start=tic;
        % Run optimization algorithm
        [MSE_values(i),MSE_list] = FindOptA(h, N, numNodes, sigma_sqr, rho, Pmax, L,O,eps);
        time(i)=toc(t_start);
    end
    mean_MSE(idx)=mean(MSE_values);
    mean_time(idx)=mean(time);
end
% Display average MSE
%disp(['Average MSE: ', num2str(mean(MSE_values))]);
figure;
plot(N_values, mean_time, '-o', 'LineWidth', 2, 'MarkerSize', 8);
%set(gca, 'YScale', 'log'); % Set y-axis to log scale
%ylim([0, 6]); % Set y-axis limits

xlabel('Number of Selected Antennas (O)');
ylabel('Mean MSE');
title('Mean MSE vs. Number of Selected Antennas');
grid on;
function [min_mse,MSE_list] = FindOptA(hk, N, K, sigma_sqr, rho, Pmax, L,O,eps) 

    a_hat = (randn(N,1) + 1i*randn(N,1)) / sqrt(2);
    A_hat = a_hat * a_hat';
    
    % Precompute Hk matrices
    for k = 1:K
        Hk(:,:,k) = hk(:,k) * hk(:,k)';
    end
    
    % Iterative optimization using DC approach
    for t=1:N-O
        % Compute subgradient
        [U, ~, ~] = svd(A_hat);
        v = U(:,1);
        subgrad = v * v';
        [~, sorted_idx] = sort(diag(A_hat), 'ascend');
        min_idx = sorted_idx(1:t);
        
        % Solve the convex subproblem for ZF (all signals aligned)
        cvx_begin quiet
            variable A(N,N) semidefinite hermitian
            minimize(sigma_sqr*trace(A) + rho*trace(A*(eye(N)-subgrad)'))
            
            subject to
                sum(diag(A(min_idx, min_idx))) <= eps
                for k = 1:K
                    real(trace(A * Hk(:,:,k))) >= 1 / Pmax;
                end
                %trace(A)>=1

        cvx_end 
        %disp(cvx_optval)
        %[U,S,~]=svd(A);
        %U1 = U(:,1);
        %sigma1 = S(1,1);
        %ahat_t = U1 * sqrt(sigma1);
        A_hat=A;
    end
    [U, Sigma] = svd(A);
    a = zeros(N, L);
    MSE = zeros(1, L);
    r = size(Sigma,1);
    for j = 1:L
        % Generate random vector
        n = (randn(r,1) + 1i*randn(r,1))/sqrt(2);
        a(:,j) = U * sqrt(Sigma) * n;
        [~, idx] = sort(abs(a(:,j)), 'ascend');  % Sort by magnitude
        a(idx(1:N-O), j) = 0;
        %[~, min_idx] = mink(a(:,j), N-O);
        % Normalize to satisfy ZF constraints
        scaling_factor = min(abs(a(:,j)'*hk)) * sqrt(Pmax);
        a(:,j) = a(:,j) / scaling_factor;
        %end
        MSE(j) = sigma_sqr * norm(a(:,j))^2;
    end

    % Find best result
    [min_mse,idx] = min(MSE);
    %ahat_t = a(:,idx);
    %A_hat=ahat_t*ahat_t';
    MSE_list=MSE;
    %end
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
