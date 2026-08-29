clear;
%% Simulation Parameters
Pmax = 10;
sigma_sqr = 1;
numTrials = 1000;       
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
N_values = 5:5:30;
mean_MSE = zeros(1, length(N_values));
for idx = 1:length(N_values)
    N = N_values(idx);
    MSE_values = zeros(1, numTrials);
    for i = 1:numTrials
        % Generate channel coefficients
        h = generateChannelCoefficients(N, numNodes, areaSize, sinkPos, ref_dist, ref_PL_dB, rx_gain_dB);
        
        % Run optimization algorithm
        MSE_values(i) = AntSelect(h, N, numNodes, sigma_sqr,Pmax);
    end
    mean_MSE(idx) = mean(MSE_values);
end
figure;
plot(N_values, mean_MSE, '-o', 'LineWidth', 2, 'MarkerSize', 8);
%set(gca, 'YScale', 'log'); % Set y-axis to log scale
ylim([0, 10]); % Set y-axis limits

xlabel('Number of Antennas (N)');
ylabel('Mean MSE [log scale]');
title('Mean MSE vs. Number of Antennas');
grid on;
function [min_MSE,idx] = AntSelect(hk, N, K, sigma_sqr,Pmax)
    MSE = zeros(1, N);  
    g = zeros(N, K);
    for i=1:N
        h_sort = sort(abs(hk(i,:)), 'ascend');
        h_sum = cumsum(h_sort);
        h_sumsqr = cumsum(h_sort.^2);
        g(i,:) = sqrt(Pmax) * h_sum ./ (sigma_sqr + Pmax * h_sumsqr);
        [gi,idx]=max(g(i,:));
        if gi>1/(h_sort(idx)*sqrt(Pmax))
            a=1/(h_sort(idx)*sqrt(Pmax));
        elseif gi<1/(h_sort(idx+1)*sqrt(Pmax)) && idx < K
            a=1/(h_sort(idx+1)*sqrt(Pmax));
        else 
            a=gi;
        end

        for k = 1:idx
            MSE(i) = MSE(i) + abs(a*h_sort(k)*sqrt(Pmax)-1)^2;
        end
        MSE(i) = MSE(i) + sigma_sqr*a^2;
    end
    
    [min_MSE, idx] = min(MSE);
end
%{
function h = generateChannelCoefficients(N, numNodes, areaSize, sinkPos, freq, antennaHeight, ref_dist, ref_PL_dB, rx_gain_dB)
    % Constants
    c = 3e8;  % Speed of light (m/s)
    lambda = c / freq;  % Wavelength (m)
    
    % Generate random node positions in a square area
    nodes = rand(numNodes, 2) * areaSize;
    
    % Compute Euclidean distances to the sink
    distances = sqrt((nodes(:,1) - sinkPos(1)).^2 + (nodes(:,2) - sinkPos(2)).^2);
    distances = max(distances, lambda/(4*pi)); % Prevent division by zero or very small values
    
    % Calculate crossover distance for hybrid model
    d_cross=90;
    PL_dB = zeros(1, numNodes);
    
    % Compute path loss based on distance range
    is_free_space = distances <= d_cross;
    
    % Free-space path loss model
    PL_dB(is_free_space) = 20 * log10(distances(is_free_space)) + 20 * log10(4 * pi / lambda);
    
    % Two-ray ground reflection model
    PL_dB(~is_free_space) = 40 * log10(distances(~is_free_space)) - 20 * log10(antennaHeight^2) ;
    
    % Compute uncalibrated reference path loss at ref_dist
    if ref_dist <= d_cross
        ref_PL_uncalibrated = 20 * log10(ref_dist) + 20 * log10(4 * pi / lambda);
    else
        ref_PL_uncalibrated = 40 * log10(ref_dist) - 20 * log10(antennaHeight^2);
    end
    
    % Apply calibration offset
    calibration_offset = ref_PL_dB - ref_PL_uncalibrated;
    PL_dB = PL_dB + calibration_offset;
    
    % Adjust for receiver gain
    effective_PL_dB = PL_dB - rx_gain_dB;
    
    % Convert path loss (dB) to linear channel gain
    channel_gain = 10.^(-effective_PL_dB / 10);
    
    % Apply Rayleigh fading
    h = (randn(N, numNodes) + 1i * randn(N, numNodes)) / sqrt(2); % Complex Gaussian
    h = h .* sqrt(reshape(channel_gain, [1, numNodes])); % Apply channel gain
    
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