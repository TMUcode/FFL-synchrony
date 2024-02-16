%THE FOLLOWING CODE COMPUTES THE SYNC. MEASURE SIGMA FOR THE L-O SYSYEM
% Clear all existing variables and close all figures
clear all; close all;

% FUNCTION PARAMETERS
% Time parameters
dt = 0.01;           % Time step
T = 200;             % Duration of simulation
t0 = 50/dt;          % Index of transient time
t = (0:dt:T);        % Time domain

% Coupling parameters
d = 0.1;             % Common coupling strength
d21 = d;             % Coupling between Osc. 1 & 2
d31 = d;             % Coupling between Osc. 1 & 3
d32 = d;             % Coupling between Osc. 2 & 3

num_trials = 50;     % Number of trials
del_vec = logspace(-4, 0.5, 20);  % Vector of noise intensities

sigma_temp = zeros(num_trials, length(del_vec)); % Temporary storage for synchronization measure

% Loop over different noise intensities
for j = 1:length(del_vec)
    % Current noise intensity
    del_1 = del_vec(j);  % Driving noise intensity for Osc. 1
    del_2 = 0.01;        % Fixed noise intensity for Osc. 2
    del_3 = 0.01;        % Fixed noise intensity for Osc. 3
    
    % Loop over trials
    for i = 1:num_trials
        % Simulate model
        [X] = L_O_Motif_simulation(del_1, del_2, del_3, d21, d31, d32, dt, T);

        % Smooth time series (after removing transient time)
        x1 = smoothdata(X(1, t0:end), 'gaussian', 100);
        x2 = smoothdata(X(3, t0:end), 'gaussian', 100);
        x3 = smoothdata(X(5, t0:end), 'gaussian', 100);

        % Compute average peak amplitude for each oscillator
        [pks1, ~] = findpeaks(x1); amp_1 = mean(pks1);
        [pks2, ~] = findpeaks(x2); amp_2 = mean(pks2);
        [pks3, ~] = findpeaks(x3); amp_3 = mean(pks3);

        % Normalize each oscillator's time series by its average peak amplitude
        X_x = [x1/amp_1; x2/amp_2; x3/amp_3];

        % Calculate standard deviation of the normalized time series for each neuron
        sig_neuron = std(X_x);

        % Average standard deviation across neurons
        sigma_temp(i, j) = mean(sig_neuron);
    end % end trials
end % end noise intensities

% Compute mean synchronization measure across trials
sigma = mean(sigma_temp);

% Plot synchronization measure against noise intensity
figure(1)
semilogx(del_vec, sigma, '*-k')
xlabel('\delta_1')
ylabel('\sigma')
title('Network Synchronization')

% Compute and display mean synchronization measure
mean_sigma = mean(sigma);
disp(['Mean synchronization measure: ', num2str(mean_sigma)]);
