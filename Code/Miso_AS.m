clear;
%% Simulation Parameters
N = 10;
O=5;
eps=0.001;
Pmax = 10;
sigma_sqr = 1;
rho = 2;
beta=2;
numTrials = 100;       
numNodes = 50;          
sinkPos = [100, 0];
areaSize = 200;             
freq = 2.4e9;
c = 3e8;
L = 1000;
ref_dist = 90;           
ref_PL_dB = 80;          
rx_gain_dB = 70;         

%% Simulation
N_values = 10;
mean_MSE = zeros(1, length(N_values));
for idx = 1:length(N_values)
    disp(idx)
    N = N_values(idx);
    MSE_values = zeros(1, numTrials);
    for i = 1:numTrials
        % Generate channel coefficients
        h = generateChannelCoefficients(N, numNodes, areaSize, sinkPos, ref_dist, ref_PL_dB, rx_gain_dB);
        
        % Run optimization algorithm
        MSE_values(i) = miso_as(h, N, numNodes, sigma_sqr, rho, Pmax, L,O,eps,beta);
    end
    mean_MSE(idx) = mean(MSE_values);
end
figure;
plot(N_values, mean_MSE, '-o', 'LineWidth', 2, 'MarkerSize', 8);
%ylim([6, 17]); % Set y-axis limits

xlabel('Number of Antennas (N)');
ylabel('Mean MSE [log scale]');
title('Mean MSE vs. Number of Antennas');

function min_mse = miso_as(hk, N, K, sigma_sqr, rho, Pmax, L,O,zi,beta) 

    a_hat = (randn(N,1) + 1i*randn(N,1)) / sqrt(2);
    A_t_prev = a_hat * a_hat';
    
    % Precompute Hk matrices
    for k = 1:K
        Hk(:,:,k) = hk(:,k) * hk(:,k)';
    end
    num_Si=0;

    % Iterative optimization using DC approach
    for t=1:N-O

        scores=abs(a_hat'*hk);
        [~,top_sorted]=sort(scores,'ascend');
        top_sorted=top_sorted(1:num_Si);
        not_S = setdiff(1:K, top_sorted);
        % Compute subgradient
        [U, ~, ~] = svd(A_t_prev);
        v = U(:,1);
        subgrad = v * v';

        [~, sorted_idx] = sort(diag(A_t_prev), 'ascend');
        min_idx = sorted_idx(1:t);
        
        % Solve the convex subproblem for ZF (all signals aligned)
        %cvx_solver sedumi
        %eps = zi * trace(A_t_prev);
        cvx_begin quiet
            variable A(N,N) semidefinite
            miso_tr = 0;
            if ~isempty(top_sorted)
                for i = 1:length(top_sorted)
                    miso_tr = miso_tr + real(trace(A * Hk(:,:,top_sorted(i))));
                end
            end
            minimize(sigma_sqr * trace(A) + rho * trace(A * (eye(N) - subgrad))-beta * miso_tr) 
            subject to
                % Constraints for k in S
                %if ~isempty(top_sorted)
                %    for k = 1:length(top_sorted)
                %        real(trace(A * Hk(:,:,top_sorted(k)))) <= 1 / Pmax;
                %    end
                %end
                % Constraints for k not in S
                if ~isempty(not_S)
                    for k = 1:length(not_S)
                        real(trace(A * Hk(:,:,not_S(k)))) >= 1 / Pmax;
                    end
                end
                sum(diag(A(min_idx, min_idx))) <= eps
        cvx_end 
        if strcmp(cvx_status, 'Solved')
            disp("Great")
            disp(cvx_optval)
        else
                disp(cvx_status)
                A = A_t_prev;
        end
        [U, Sigma,~] = svd(A);
        a = zeros(N, L);
        MSE = zeros(1, L);
        r = size(Sigma,1);
        b=zeros(L,K);
        for j = 1:L
            % Generate random vector
            n = (randn(r,1) + 1i*randn(r,1))/sqrt(2);
            a(:,j) = U * sqrt(Sigma) * n;

            [~, idx] = sort(abs(a(:,j)), 'ascend');  % Sort by magnitude
            %if t > 1
            a(idx(1:min(t, N-O)), j) = 0;
            %end
            remaining_nodes = setdiff(1:K, top_sorted);

            if ~isempty(remaining_nodes)
                min_val = min(abs(a(:,j)' * hk(:, remaining_nodes)));
            else
                min_val = 1; % Default nonzero value to avoid division by zero
            end
            a(:,j) = (a(:,j) / min_val) * (1 / sqrt(Pmax));
            % Compute b_k for each node to achieve alignment
            for k = 1:length(top_sorted)
                MSE(j)=MSE(j)+(1-abs(a(:,j)' * hk(:,top_sorted(k))*sqrt(Pmax)))^2;
            end
            
            % For ZF, MSE is just the noise term since all signals are aligned
            MSE(j) = MSE(j)+sigma_sqr * norm(a(:,j))^2;
        end
        
        [min_mse,min_idx] = min(MSE);
        a_hat=a(:,min_idx);
        A_t_prev=A;
        num_Si=num_Si+1;
    end
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
