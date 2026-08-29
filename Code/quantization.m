clear;
%% Simulation Parameters
N = 10;
O=1;
eth = 0.001;
eps=0.001;
Pmax = 10;
sigma_sqr = 1;
rho = 2;
numTrials = 10;       
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
Nq_values = 5:5:30;
mean_MSE = zeros(1, length(Nq_values));
mean_MSE_th = zeros(1, length(Nq_values));
mean_MSE_opt = zeros(1, length(Nq_values));
for idx = 1:length(Nq_values)
    MSE_values = zeros(1, numTrials);
    MSE_theory = zeros(1, numTrials);
    MSE_opt = zeros(1, numTrials);
    for i = 1:numTrials
        disp(i)
        % Generate channel coefficients
        h = generateChannelCoefficients(N, numNodes, areaSize, sinkPos, ref_dist, ref_PL_dB, rx_gain_dB);
        [MSE_values(i),MSE_theory(i),MSE_opt(i)] = ZF(h, N, numNodes, sigma_sqr, rho, eth,Pmax, L,O,Nq_values(idx));
    end
    mean_MSE(idx)=mean(MSE_values);
    mean_MSE_th(idx)=mean(MSE_theory);
    mean_MSE_opt(idx)=mean(MSE_opt);
end
% Display average MSE
%disp(['Average MSE: ', num2str(mean(MSE_values))]);
figure;
plot(Nq_values, mean_MSE, '-o', 'LineWidth', 2, 'MarkerSize', 8);hold on;
plot(Nq_values, mean_MSE_th, '-^', 'LineWidth', 2, 'MarkerSize', 8);
%plot(Nq_values, mean_MSE_opt, '-v', 'LineWidth', 2, 'MarkerSize', 8);hold off;
%set(gca, 'YScale', 'log'); % Set y-axis to log scale
%ylim([0, 6]); % Set y-axis limits

xlabel('Quantization Number (Nq)');
ylabel('Mean MSE');
title('Mean MSE vs. Quantization Number');
grid on;

function [MSE_ZF,MSE_th,MSE_no_q] = ZF(hk, N, K, sigma_sqr, rho, eth, Pmax, L,O,Q) 
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
    eps=exp(1j*angle(a_min))*(norm(a_min,1)/N)-a_min;
    MSE_no_q=sigma_sqr*norm(eps+a_min)^2;
    Q_mat=zeros(N,N);
    for k = 1:K
        inner_prod = a_min' * hk(:,k); % Compute a^H h_k
        b(k) = 1 / inner_prod;
        Q_mat=Q_mat+abs(b(k))^2*(hk(:,k)*hk(:,k)');
        MSE_no_q=MSE_no_q+abs(eps'*hk(:,k)*b(k))^2;
    end
    [~,~,a_hybrid]=unit_modulus_factorization(a_min,O,Q);
    
    %MSE_ZF=compute_mse(a_hybrid,b,hk,K,sigma_sqr);
    %lambda_max = eigs(Q_mat, 1, 'largestreal');
    vals = eigs(Q_mat,N, 'largestreal');  % returns top O+1 eigenvalues (unsorted)
    vals = sort(vals, 'descend');
    %if N~=O
    %    lambda_max = vals(O+1);
    %    MSE_th = (norm(a_min)^2-((norm(a_min,1)*cos(pi/Q))^2)/N) * (real(lambda_max) - sigma_sqr) + sigma_sqr*norm(a_min)^2;
    %else
    %    MSE_th=sigma_sqr*norm(a_min)^2;
    %end
    MSE_th=norm(a_min)^2-norm(a_min,1)^2*(cos(pi/Q))^2/N;
    MSE_ZF=norm(a_min-a_hybrid)^2;
end


function MSE = compute_mse(a, b, hk, K, sigma_sqr)
    MSE = 0;
    for k = 1:K
        MSE = MSE + abs(a' * hk(:,k) * b(k)-1)^2;
    end
    MSE = MSE+sigma_sqr * norm(a)^2;
end
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
