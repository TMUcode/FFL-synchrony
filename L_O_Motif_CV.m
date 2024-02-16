clear all; close all;
% CODE DESCRIPTION:
% This code computes the coefficient of variation (CV) of the inter-spike
% intervals (ISIs) for oscillator 3 in a motif of coupled oscillators.

% time
dt = 0.01;       % Time step
T = 200;         % Duration of simulation
t0 = 50/dt;      % Index of transient
t = (0:dt:T);    % Time domain
t1 = t(t0:end);  % Time domain after transient period

% Coupling strengths
d21 = 0.1;  % Coupling between Oscillator 1 & 2
d31 = 0.1;  % Coupling between Oscillator 1 & 3
d32 = 0.1;  % Coupling between Oscillator 2 & 3

num_trials = 100;   % Number of trials
del_vec = logspace(-3,0.5,20); % Vector of noise intensities


for j = 1: length(del_vec) % Loop over different noise intensities


    % Noise intensity
    del_1 = del_vec(j);  % Driving noise intensity for Oscillator 1
    del_2 = 0.01;        % Fixed noise intensity for Oscillator 2
    del_3 = 0.01;        % Fixed noise intensity for Oscillator 3


    % Initialize variables to store spike times and ISIs for Oscillator 3
    ST3 = [];   % Spike times
    ISI3 = [];  % Inter-spike intervals

    for i = 1:num_trials  % Perform multiple trials for each noise intensity

        % SIMULATE MODEL
        [X] = L_O_Motif_simulation(del_1,del_2,del_3,d21,d31,d32,dt,T);

        % Extract data for Oscillator 3 and smooth the time series
        x3 = smoothdata(X(5, t0:end), 'gaussian', 100);

        %Find peaks smooted data (un-normzalized)
        [pks3,~] = findpeaks(x3); amp3 = mean(pks3);

        x3_norm = x3/amp3; % Amplitude normalized timeseries

        % Find peaks in the smoothed data to identify spike times (normalized)
        [pks3, locs3] = findpeaks(x3_norm, 'MinPeakDistance', 3,  ...
            'MinPeakHeight', 0, 'MinPeakProminence', 0.05);

        % Calculate spike times and inter-spike intervals for trial i
        st3 = t1(locs3);           % Spike times for current trial
        isi3 = diff(st3);          % Inter-spike intervals for trial i
        % Collect inter-spike intervals across trials
        ISI3 = [ISI3, isi3];
    end

    % Calculate coefficient of variation (CV) for ISI distribution of Oscillator 3
    r3(j) = std(ISI3) / mean(ISI3);

end

% Plot CV of ISI distribution against noise intensity
semilogx(del_vec, r3, ':ok')
xlabel('\delta_1')
ylabel('CV')
title('Coefficient of Variation (CV) of Inter-Spike Intervals (ISIs) for Oscillator 3')