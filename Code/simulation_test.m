clear;
%% Simulation Parameters
Pmax = 10;
sigma_sqr = 1;
numTrials = 1000;       
numNodes = 50;
N_values = 20:20:200;
O=4;
N=50;
Q=10;
%% Simulation
mean_MSE_sim = zeros(1, length(N_values));
mean_MSE_th = zeros(1, length(N_values));
mean_MSE_opt = zeros(1, length(N_values));
for idx = 1:length(N_values)
    MSE_sim = zeros(1, numTrials);
    MSE_theory = zeros(1, numTrials);
    MSE_opt = zeros(1, numTrials);
    numNodes=N_values(idx);
    %numNodes=N_values(idx);
    for i = 1:numTrials
        if mod(i, 200) == 0
            fprintf('i = %d, idx = %d\n', i, idx);
        end
        % Generate channel coefficients
        h = generateChannelCoefficients(N, numNodes);
        [MSE_sim(i),MSE_theory(i),MSE_opt(i)] = sim(h, N, numNodes, sigma_sqr,Pmax,O,Q);
    end
    mean_MSE_sim(idx)=mean(MSE_sim);
    mean_MSE_th(idx)=mean(MSE_theory);
    mean_MSE_opt(idx)=mean(MSE_opt);
end

figure;
set(gcf, 'Units', 'inches', 'Position', [1, 1, 6, 3.5]);
%plot(N_values, mean_MSE_sim, '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
plot(N_values, mean_MSE_th, '-^', 'LineWidth', 2, 'MarkerSize', 8);hold on;
plot(N_values, mean_MSE_opt, '-v', 'LineWidth', 2, 'MarkerSize', 8);hold off;
xticks([40 80 120 160 200]);
xlabel('Number of Nodes (K)');
ylabel('Mean MSE');
legend('MSE upper bound','Ideal MSE', 'Location', 'northwest');
%title('Mean MSE vs.  Number of Nodes');

%grid on;
set(gca, 'FontSize', 12);  % axes ticks font size
grid off;

% Save as EPS
print(gcf, 'mse_bound', '-depsc2');
function [MSE_sim,MSE_upper,MSE_opt]=sim(h,N,K,sigma_sqr,Pmax,O,Q)
    inverse_term=sigma_sqr*eye(N);
    for k=1:K
        inverse_term=inverse_term+Pmax*h(:,k)*h(:,k)';
    end
    second_term=zeros(N,1);
    b=zeros(K,1);
    for k = 1:K
        theta = 2 * pi * rand();                          
        b(k) = sqrt(Pmax) * exp(1j * theta);              
        second_term = second_term + h(:,k) * b(k);  
    end
    a_opt=inverse_term\second_term;
    a_hat=unit_modulus_factorization(a_opt,O,Q);
    MSE_opt=compute_mse(a_opt,b,h,K,sigma_sqr);
    %MSE_sim=norm(a_opt-a_hat)^2;
    MSE_upper=(K*Pmax+sigma_sqr)*N*(2*sin(pi/(2*Q)))^O*max(abs(a_opt))^2+MSE_opt;
    MSE_sim=compute_mse(a_hat,b,h,K,sigma_sqr);
    %MSE_th=(K*Pmax+sigma_sqr)*norm(a_hat-a_opt)^2+MSE_opt;
    %bound=N*(2*sin(pi/(2*Q)))^O*max(abs(a_opt))^2;
    %MSE_th=(K*Pmax+sigma_sqr)*bound+MSE_opt;
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
        
        [~, idx1] = min(abs(w1 - W_values));
        [~, idx2] = min(abs(w2 - W_values));
        
        W(i,1) = W_values(idx1);
        W(i,2) = W_values(idx2);
    end 
end

function h = generateChannelCoefficients(N, numNodes)
    h = (randn(N, numNodes) + 1i*randn(N, numNodes))/sqrt(2);  % CN(0,1)
end