clear all; close all;

% CODE DESCRIPTION: The following code computes the spike times and ISIs
% of the oscillators x_i

% FUNCTION PARAMETERS
% time
dt = 0.01;           % timestep
T = 200;             % duration of simulation
t0 = 50/dt;          % transient time (index of transient)
t = (0:dt:T);        % time domain
t1 = t(t0:end);

% coupling
d21 = 0.1;           % coupling between Osc. 1 & 2
d31 = 0.1;           % coupling between Osc. 1 & 3
d32 = 0.1;           % coupling between Osc. 2 & 3

num_trials = 25;     % number of trials
del_vec = logspace(-3,0.5,20); % vector of noise intensities

% Loop over different noise intensities
for j = 1:length(del_vec)

    % Noise intensity for each oscillator
    del_1 = del_vec(j);  
    del_2 = del_vec(j);  
    del_3 = del_vec(j);  

    % Perform multiple trials for each noise intensity
    for i = 1:num_trials

        % SIMULATE MODEL
        [X] = L_O_Motif_simulation(del_1, del_2, del_3, d21, d31, d32, dt, T);

        % Smooth the data
        x1 = smoothdata(X(1, t0:end), 'gaussian', 100);
        y1 = smoothdata(X(2, t0:end), 'gaussian', 100);
        x2 = smoothdata(X(3, t0:end), 'gaussian', 100);
        y2 = smoothdata(X(4, t0:end), 'gaussian', 100);
        x3 = smoothdata(X(5, t0:end), 'gaussian', 100);
        y3 = smoothdata(X(6, t0:end), 'gaussian', 100);

        % Find peaks and calculate their average amplitude
        [pks1, ~] = findpeaks(x1, 'MinPeakDistance', 3, 'MinPeakHeight', 0); 
        amp1(i) = mean(pks1); 
        [pks2, ~] = findpeaks(x2, 'MinPeakDistance', 3, 'MinPeakHeight', 0); 
        amp2(i) = mean(pks2);
        [pks3, ~] = findpeaks(x3, 'MinPeakDistance', 3, 'MinPeakHeight', 0); 
        amp3(i) = mean(pks3);

        % Calculate the average magnitude of oscillation
        r1 = sqrt(x1.^2 + y1.^2); 
        r1_avg(i) = mean(r1);
        r2 = sqrt(x2.^2 + y2.^2); 
        r2_avg(i) = mean(r2);
        r3 = sqrt(x3.^2 + y3.^2);  
        r3_avg(i) = mean(r3); 

    end
    
    % Calculate mean amplitudes and mean magnitudes for each noise intensity
    A1(j) = mean(amp1);
    A2(j) = mean(amp2);
    A3(j) = mean(amp3);

    R1(j) = mean(r1_avg);
    R2(j) = mean(r2_avg);
    R3(j) = mean(r3_avg);

end

% Plot the mean magnitudes against noise intensity for each oscillator
semilogx(del_vec, R1, '--k', 'LineWidth', 1.5)
hold on
semilogx(del_vec, R2, ':k', 'LineWidth', 1.5)
semilogx(del_vec, R3, 'k', 'LineWidth', 1.5)
legend('x_1', 'x_2', 'x_3')
