clc; close all; clear all;

load('data1.mat');                 % Import samples' file
whos                           % Check the contents of the loaded .mat file

% Plot the signals
figure(1);
plot(y, 'black', 'LineWidth', 1.3); 
hold on; grid on; 
plot(u, 'r', 'LineWidth', 1.7); 
hold off; 

% Add title and labels
title('Signals y[n] and u[n]');
xlabel('Sample Index');
ylabel('Signal Values');

% Add legend
legend('Output Signal y[n]', 'Input Signal u[n]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
orders = [2, 3, 4];  % Orders of the unknown filter

% Define specific mu values for each order
% Step size values for the LMS algorithm
mu_values_all = {
    [0.05, 0.01],   % For p = 2
    [0.15, 0.05],   % For p = 3
    [0.3, 0.03]     % For p = 4
};

% Calculate the sum of squared values of input signal u
N = length(u);
sum_sq_u = sum(abs(u).^2);

% Define colors and line styles
colors = {'#0EC3EB', '#FA288C'};  % Colors for different mu values
line_styles = {'-', '--', ':', '-.'};  % Line styles for different coefficients

fprintf('The step size values must meet the condition ...\n');

% Create a figure to plot all J(w) values in subplots
figure_MSE = figure;
subplot_idx = 1;

% Execute the adaptive LMS algorithm for each combination of order p and step size mu
for p_idx = 1:length(orders)
    p = orders(p_idx);
    mu_values = mu_values_all{p_idx};
    
    % Compute upper bound for mu based on the given formula
    mu_upper_bound = 2 * N / ((p + 1) * sum_sq_u);
    
    % Print mu values bound
    fprintf('For order p = %d:\n', p);
    fprintf('0 < μ < %.4f\n', mu_upper_bound);

    % Initialize array to store J values
    J_values = zeros(length(mu_values), length(y) - p + 1);
    
    % Create a new figure for each order
    figure; hold on; grid on;
    
    for mu_idx = 1:length(mu_values)
        mu = mu_values(mu_idx);

        % Initialize the adaptive filter coefficients
        W = zeros(p, 1);  

        % Initialize a matrix to store coefficients at each iteration
        W_all = zeros(p, length(y)-p+1);
        
        % Initialize the error array to compute J
        e_all = zeros(1, length(y)-p+1);
        
        % Apply the LMS algorithm
        for n = p:length(y)
            % Filter input (signal u)
            u_input = u(n:-1:n-p+1).';  % Ensure u_input is a column vector
            
            % Compute the output of the adaptive filter
            y_hat = W' * u_input;  % Matrix multiplication to calculate the estimated output
            
            % Compute the error (difference between the signal y[n] and the filter output)
            e = y(n) - y_hat;
            
            % Store the error
            e_all(n-p+1) = e;
            
            % Update the filter coefficients based on the error and the input
            W = W + mu .* e .* conj(u_input);
            
            % Store the current coefficients
            W_all(:, n-p+1) = W;
        end
        
        % Compute the mean squared error J for each iteration and store it
        J_values(mu_idx, :) = abs(e_all).^2;
        
        % Plot the convergence of the adaptive filter coefficients
        for coeff_idx = 1:p
            plot(W_all(coeff_idx, :), 'Color', colors{mu_idx}, 'LineStyle', line_styles{coeff_idx}, ...
                'LineWidth', 2, 'DisplayName', sprintf('w_n[%d] for \\mu=%.2f', coeff_idx-1, mu));
        end
    end
    
    % Add title and labels for each plot
    title(sprintf('Adaptive Filtering via LMS with p=%d', p));
    xlabel('Iteration n');
    ylabel('LMS-Based FIR Filter Coefficients');
    legend show;
    ylim([-3 4]);  
    hold off;
    
    % Convert J_values to dB
    J_values_dB = 10 * max(log10(J_values), 0); 
    
    % Plot J values in dB for each mu value in subplots
    figure(figure_MSE);
    subplot(3, 1, subplot_idx);  % Create subplots for each order
    hold on; grid on;
    for mu_idx = 1:length(mu_values)
        plot(1:length(J_values_dB(mu_idx, :)), J_values_dB(mu_idx, :), 'Color', colors{mu_idx}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('p=%d, \\mu=%.2f', p, mu_values(mu_idx)));
    end
    title(sprintf('MSE Convergence of J(w) for p=%d', p));
    xlabel('Iteration n');
    ylabel('MSE in dB');
    legend('Location', 'best');
    hold off;
    
    subplot_idx = subplot_idx + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End of file                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
