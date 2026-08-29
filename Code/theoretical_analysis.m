clear;
%% Simulation Parameters
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
eps=10^-2;
N_values = 2:2:10;
O=4;
N=20;
Q=10;
%% Simulation
mean_MSE = zeros(1, length(N_values));
mean_MSE_th = zeros(1, length(N_values));
mean_MSE_opt = zeros(1, length(N_values));
mean_time=zeros(1,length(N_values));
for idx = 1:length(N_values)
    MSE_values = zeros(1, numTrials);
    MSE_theory = zeros(1, numTrials);
    MSE_opt = zeros(1, numTrials);
    time=zeros(1,numTrials);
    O=N_values(idx);
    %numNodes=N_values(idx);
    for i = 1:numTrials
        if mod(i, 20) == 0
            fprintf('i = %d, idx = %d\n', i, idx);
        end
        % Generate channel coefficients
        h = generateChannelCoefficients(N, numNodes, areaSize, sinkPos, ref_dist, ref_PL_dB, rx_gain_dB);
        t_start=tic;
        [MSE_values(i),MSE_opt(i),MSE_theory(i)] = ZF(h, N, numNodes, sigma_sqr, rho, eth,Pmax, L,O,Q);
        time(i)=toc(t_start);
    end
    mean_MSE(idx)=mean(MSE_values);
    mean_MSE_th(idx)=mean(MSE_theory);
    mean_MSE_opt(idx)=mean(MSE_opt);
    mean_time(idx)=mean(time);
end
% Display average MSE
%disp(['Average MSE: ', num2str(mean(MSE_values))]);
figure;
%semilogy(N_values, mean_MSE, '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
%emilogy(N_values, mean_MSE_th, '-^', 'LineWidth', 2, 'MarkerSize', 8);
plot(N_values, mean_time); 
%plot(N_values, mean_MSE, '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
%plot(N_values, mean_MSE_th, '-^', 'LineWidth', 2, 'MarkerSize', 8);
%plot(N_values, mean_MSE_opt, '-v', 'LineWidth', 2, 'MarkerSize', 8);hold off;

xlabel('Number of Antennas (N)');
ylabel('Mean MSE');
title('Mean MSE vs. Number of Antennas');
grid on;

function [MSE_ZF,MSE_opt,MSE_th] = ZF(hk, N, K, sigma_sqr, rho, eth, Pmax, L,O,Q) 
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
    a_hybrid=unit_modulus_factorization(a_min,O,Q);
    eps=a_hybrid-a_min;
    %MSE_th=sigma_sqr*norm(a_min)^2;
    delta = (2 * sin(pi / (2 * Q)))^(O / 2);
    aopt_max = max(abs(a_min));
    aopt_min=min(abs(a_min));
    for k = 1:K
        inner_prod = a_min' * hk(:,k); % Compute a^H h_k
        b(k) = 1 / inner_prod;
        %MSE_th=MSE_th+abs(b(k))^2*norm(eps)^2;

    end
    MSE_th=delta^2*norm(a_min)^2;
    MSE_ZF=norm(eps)^2;
    %MSE_ZF=compute_mse(a_hybrid,b,hk,K,sigma_sqr);
end


function MSE = compute_mse(a, b, hk, K, sigma_sqr)
    MSE = 0;
    for k = 1:K
        MSE = MSE + abs(a' * hk(:,k) * b(k)-1)^2;
    end
    MSE = MSE+sigma_sqr * norm(a)^2;
end

function an = unit_modulus_factorization(v, O, Q)
    an=zeros(length(v),1);
    residue=v;
    t=floor(O/2);
    for i=1:t
        [ao,W]=approx(residue,Q);
        an=an+W*ao;
        residue=residue-W*ao;
    end
    theta_Q = 2*pi*(0:Q-1)/Q;
    W_values = exp(1i * theta_Q);  
    N = length(residue);
    w_best = zeros(N, 1);
    if mod(O,2) ~=0
        for i = 1:N
            phase_i = angle(residue(i));
            [~, idx] = min(abs(wrapToPi(theta_Q - phase_i)));  % Nearest quantized phase
            w_best(i) = W_values(idx);
        end

        % Compute optimal ao for this w
        ao = (w_best' * residue) / N;
        residue = residue - ao * w_best;
        an=an+ao*w_best;
    end
end
function [u,W]=approx(v,Q)
    theta_Q = 2*pi*(0:(Q-1))/Q;
    W_values = exp(1i * theta_Q);
    N=length(v);
    u=zeros(2,1);
    u(1)=(max(abs(v))+min(abs(v)))/2;
    u(2)=(max(abs(v))-min(abs(v)))/2;
    W=zeros(N,2);
    for i=1:N
        a1 = abs(v(i))^2 + max(abs(v)) * min(abs(v));
        b1 = abs(v(i)) * (max(abs(v)) + min(abs(v)));
        theta_i = angle(v(i)) - acos(min(max(a1 / b1, -1), 1));  % Clamped inside [-1,1] to avoid domain error
        w1 = exp(1j * theta_i);

        a2 = abs(v(i))^2 - max(abs(v)) * min(abs(v));
        b2 = abs(v(i)) * (max(abs(v)) - min(abs(v)));
        theta_i = angle(v(i)) + acos(min(max(a2 / b2, -1), 1));  % Clamped inside [-1,1] to avoid domain error
        w2 = exp(1j * theta_i);
        idx1 = mod(round(mod(angle(w1), 2*pi) / (2*pi) * Q), Q) + 1;
        idx2 = mod(round(mod(angle(w2), 2*pi) / (2*pi) * Q), Q) + 1;
        %[~, idx1] = min(abs(w1 - W_values));
        %[~, idx2] = min(abs(w2 - W_values));
        
        W(i,1) = W_values(idx1);
        W(i,2) = W_values(idx2);
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