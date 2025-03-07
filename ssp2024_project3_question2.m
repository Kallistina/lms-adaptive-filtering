%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%               Stochastic Signal Processing 2023-2024              %
%                             Project 3                             %
%                                                                   %
%                            QUESTION 2                             %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Question 3.2
clc; close all; clear all;

load('data2.mat');                 % Import samples' file
whos                           % Check the contents of the loaded .mat file

% Plot the signals
figure(1);
plot(x, 'b', 'LineWidth', 1.3); 
hold on; grid on; 
plot(y, 'Color', "#079451");
hold off; 

% Add title, labels and legend
title('Signals x[n] and y[n]');
xlabel('Sample Index');
ylabel('Signal Values');
legend('Input Signal x[n]', 'Output Signal y[n]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define and line styles
line_styles = {'-', '--', ':', '-.'};  % Line styles for different coefficients

% Order of the unknown filter
p = 3;  

% Step size value for the LMS algorithm
mu = 0.08;

% Calculate the sum of squared values of input signal x
N = length(x);
sum_sq_x = sum(abs(x).^2);

% Compute upper bound for mu based on the given formula
mu_upper_bound = 2 * N / ((p + 1) * sum_sq_x);
fprintf('The step size values must meet the condition ...\n');

% Print mu values bound
fprintf('For order p = %d:\n', p);
fprintf('0 < μ < %.4f\n', mu_upper_bound);


% Execute the adaptive LMS algorithm for order p and step size mu

% Initialize array to store J values
J_values = zeros(1, length(y) - p + 1);
    
% Initialize the adaptive filter coefficients
W = zeros(p, 1);  

% Initialize a matrix to store coefficients at each iteration
W_all = zeros(p, length(y)-p+1);
        
% Initialize the error array to compute J
e_all = zeros(1, length(y)-p+1);

% Initialize the estimated output signal
y_hat_all = zeros(1, length(y));
        
% Apply the LMS algorithm
for n = p:length(y)
     % Filter input (signal x)
     x_input = x(n:-1:n-p+1).';  % Ensure x_input is a column vector
            
      % Compute the output of the adaptive filter
      y_hat = W' * x_input;  % Matrix multiplication to calculate the estimated output
            
      % Store the estimated output
      y_hat_all(n) = y_hat;
            
      % Compute the error (difference between the signal y[n] and the filter output)
      e = y(n) - y_hat;
            
      % Store the error
      e_all(n-p+1) = e;
            
      % Update the filter coefficients based on the error and the input
      W = W + mu .* e .* conj(x_input);
            
      % Store the current coefficients
      W_all(:, n-p+1) = W;
 end
        
% Compute the mean squared error J for each iteration and store it
J_values(1, :) = abs(e_all).^2;
figure;

% Plot the convergence of the adaptive filter coefficients
for coeff_idx = 1:p
    plot(W_all(coeff_idx, :), 'Color', '#8C71F0', 'LineStyle', line_styles{coeff_idx}, ...
        'LineWidth', 2, 'DisplayName', sprintf('w_n[%d] for \\mu=%.2f', coeff_idx-1, mu));
    hold on; grid on;
end
   
% Add title and labels for each plot
title(sprintf('Adaptive Filtering via LMS with p=%d', p));
xlabel('Iteration n');
ylabel('LMS-Based FIR Filter Coefficients');
legend show;
ylim([-3 4]);  
hold off; 
        
% Compute the mean squared error for this order and mu
mse = mean(abs(e_all).^2);
        
% Print the mean squared error for this order and mu
fprintf('MSE for order p = %d and μ = %.2f: %.4f\n', p, mu, mse);

% Plot the original and estimated output signals
figure;
plot(y, 'Color', '#F93C57', 'LineWidth', 2); 
hold on; grid on; 
plot(y_hat_all, 'Color', '#4CC1E0', 'LineWidth', 1.4);
hold off; 

% Add title and labels
title(sprintf('Signals y[n] and $\\hat{y}[n]$ for p = %d and $\\mu$ = %.2f', p, mu), 'Interpreter', 'latex');
xlabel('Sample Index');
ylabel('Signal Values');

% Add legend
legend('Output Signal y[n]', 'Estimated Output Signal y_{hat}[n]');

% Plot J values in dB
J_values_dB = 10 * max(log10(J_values), 0); 

figure;
plot(1:length(J_values_dB), J_values_dB, 'Color', "#079451", 'LineWidth', 1.5, ...
    'DisplayName', sprintf('p=%d, \\mu=%.2f', p, mu));
title(sprintf('MSE Convergence of J(w) for p=%d', p));
xlabel('Iteration n');
ylabel('MSE in dB');
legend('Location', 'best');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           End of file                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
