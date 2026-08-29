clear;
%% Simulation Parameters
N = 4;
eth = 0.001;
Pmax = 10;
sigma_sqr = 1;
rho = 1;
numTrials =10;       
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
MSE_miso = zeros(1, 2);
%% Simulation
for num_S=5:6
    MSE_values = zeros(1, numTrials);
    disp(num_S)
    for i = 1:numTrials
        % Generate channel coefficients
        h = generateChannelCoefficients(N, numNodes, areaSize, sinkPos, ref_dist, ref_PL_dB, rx_gain_dB);
        
        % Run optimization algorithm
        MSE_values(i) = FindOptA(h, N, numNodes, sigma_sqr, rho, eth, Pmax, L,num_S);
    end

    MSE_miso(num_S-4)=mean(MSE_values);
end
disp(['Min MSE: ', num2str(min(MSE_miso))]);
plot(5:6, MSE_miso, '-o', 'LineWidth', 2, 'MarkerSize', 8);
function A_ret = FindOptA(hk, N, K, sigma_sqr, rho, eth, Pmax, L,num_S) 
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

        scores=abs(a_hat'*hk);
        [~,top_sorted]=sort(scores,'ascend');
        top_sorted=top_sorted(1:num_S);
        not_S = setdiff(1:num_S, top_sorted);
        % Compute subgradient
        [U, ~, ~] = svd(A_hat);
        v = U(:,1);
        subgrad = v * v';

        % Solve the convex subproblem for ZF (all signals aligned)
        %cvx_solver sedumi
        cvx_begin quiet
            variable A(N,N) hermitian semidefinite
            expression sum_Tr_A_Hk_S
            sum_Tr_A_Hk_S = 0;
            for i = 1:length(top_sorted)
                sum_Tr_A_Hk_S=sum_Tr_A_Hk_S+real(trace(A * Hk(:,:,top_sorted(i))));
            end
            %disp(trace(A));  % Should be ~0
            %disp(trace(A * (eye(N) - subgrad)));  % Should be ~0
            %disp(imag(trace(A * Hk(:,:,1))));
            minimize(sigma_sqr * trace(A) + rho * trace(A * (eye(N) - subgrad))  +num_S - Pmax * sum_Tr_A_Hk_S)
            subject to
                % Constraints for k in S
                for k = 1:length(top_sorted)
                    real(trace(A * Hk(:,:,top_sorted(k)))) <= 1 / Pmax;
                end
                % Constraints for k not in S
                for k = 1:length(not_S)
                    real(trace(A * Hk(:,:,not_S(k)))) >= 1 / Pmax;
                end
        cvx_end
        %{
        cvx_begin quiet
            variable A(N,N) hermitian semidefinite
            expression miso_tr
            expression imag_miso
            miso_tr = -num_S/Pmax;
            imag_miso=0;
            for i = 1:length(top_sorted)
                miso_tr = miso_tr + real(trace(A * Hk(:,:,top_sorted(i))));  % Using precomputed Hk
                imag_miso = imag_miso + imag(trace(A * Hk(:,:,top_sorted(i))));
            end
            minimize(sigma_sqr * trace(A) + rho * trace(A * (eye(N) - subgrad))- Pmax * miso_tr)  
            subject to
                imag_miso==0
        cvx_end 
        %}
        % Check convergence
        disp(cvx_optval)
        MSE_curr = cvx_optval
        if abs(MSE_curr - MSE_prev)/abs(MSE_curr) <= eth 
            break
        end
        
        MSE_prev = MSE_curr;
        [U,S,~]=svd(A);
        U1 = U(:,1);
        sigma1 = S(1,1);
        ahat_t = U1 * sqrt(sigma1);
        A_hat=ahat_t*ahat_t';
    end
    
    miso_tr = 0;
    for i = 1:length(top_sorted)
        miso_tr = miso_tr + trace(A_hat*Hk(:,:,top_sorted(i)));  % Using precomputed Hk
    end
    %if(Pmax*abs(miso_tr)>num_S)
    %    A_ret=0;
    %    return;
    %end
    %disp(Pmax*miso_tr)
    % Gaussian randomization
    [U, Sigma] = svd(A);
    a = zeros(N, L);
    MSE = zeros(1, L);
    
    for j = 1:L
        % Generate random vector
        n = (randn(N,1) + 1i*randn(N,1))/sqrt(2);
        a(:,j) = U * sqrt(Sigma) * n;
        
        % Normalize to satisfy ZF constraints
        min_val = min(abs(a(:,j)' * hk(:, setdiff(1:K, top_sorted))));
        a(:,j) = (a(:,j) ./ min_val) * (1 / sqrt(Pmax));
        
        % Compute b_k for each node to achieve alignment
        for k = 1:K
            inner_prod = a(:,j)' * hk(:,k); % Compute a^H h_k
            if ismember(k, top_sorted) % If k is in S
                b(j,k) = sqrt(Pmax) * conj(inner_prod) / abs(inner_prod);
            else % If k is not in S
                b(j,k) = 1 / inner_prod;
            end
        end
        
        % For ZF, MSE is just the noise term since all signals are aligned
        %MSE(j) = sigma_sqr * norm(a(:,j))^2;
        MSE(j)=compute_mse(a(:,j),b(j,:),hk,K,sigma_sqr);
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