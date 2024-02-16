clear all; close all;
%The following script computes sync. measures for the L-O system with
%varying levels of noise and coupling

% FUNCTION PARAMETERS
% Time parameters
dt = 0.01;           % Time step
T = 200;             % Duration of simulation
t0 = 50/dt;          % Index of transient
t = (0:dt:T);        % Time domain
t1 = t(t0:end);      % Time domain after transient period

% Parameters for coupling and noise
d_vec = logspace(-4,1,20);   % Vector of coupling strengths
del_vec = logspace(-3,0,20); % Vector of noise intensities for Oscillator 1
del_2 = 0.01;                % Fixed noise intensity for Oscillator 2
del_3 = 0.01;                % Fixed noise intensity for Oscillator 3
num_trials = 25;             % Number of trials

% Loop over different coupling strengths
for k = 1:length(d_vec)
    d = d_vec(k);
    d21 = d;       % Coupling between Oscillator 1 & 2
    d31 = d;       % Coupling between Oscillator 1 & 3
    d32 = d;       % Coupling between Oscillator 2 & 3

    % Loop over different noise intensities
    for j = 1:length(del_vec)

        % Noise intensity for Oscillator 1
        del_1 = del_vec(j);

        % Initialize variables to store spike times and ISIs for Oscillator 3
        ST3 = [];
        ISI3 = [];

        % Perform multiple trials for each noise intensity
        for i = 1:num_trials

            % SIMULATE MODEL
            [X] = L_O_Motif_simulation(del_1, del_2, del_3, d21, d31, d32, dt, T);

            % Extract and smooth time series data for each oscillator
            x1 = smoothdata(X(1,t0:end),'gaussian',100);
            x2 = smoothdata(X(3,t0:end),'gaussian',100);
            x3 = smoothdata(X(5,t0:end),'gaussian',100);

            % Calculate mean amplitude for each oscillator
            [pks1,~] = findpeaks(x1); amp_1 = mean(pks1);
            [pks2,~] = findpeaks(x2); amp_2 = mean(pks2);
            [pks3,~] = findpeaks(x3); amp_3 = mean(pks3);

            % Compute sigma (standard deviation) for each neuron
            X_x = [x1/amp_1; x2/amp_2; x3/amp_3];
            sig_neuron = std(X_x);
            sigma_temp(i) = std(sig_neuron);

            % Compute R (coefficient of variation) for Oscillator 3
            x3_norm = x3/amp_3;
            [pks3,locs3] = findpeaks(x3_norm,'MinPeakDistance',3,...
                'MinPeakHeight',0,'MinPeakProminence',0.05);
            st3 = t1(locs3);     % Spike times for current trial i
            isi3 = diff(st3);    % Inter-spike intervals for trial i
            ISI3 = [ISI3, isi3]; % Cumulative ISIs

        end % End of trials loop

        % Compute mean sigma and mean R for the current noise intensity
        sigma(j) = mean(sigma_temp);
        r3(j) = std(ISI3)/mean(ISI3);

    end % End of noise intensities loop

    % Store mean sigma and mean R for each coupling strength
    sig(k,:) = sigma;
    r(k,:) = r3;
end

% Find minimum sigma and R values for each coupling strength
for coup = 1: length(d_vec)
    [min_sig_temp,loc_sig_temp] = min(sig(coup,:));
    [min_r_temp,loc_r_temp,] = min(r(coup,:));

    min_sig(coup) = min_sig_temp; loc_sig(coup) = del_vec(loc_sig_temp);
    min_r(coup) = min_r_temp; loc_r(coup) = del_vec(loc_r_temp);

end

% Plot mean square deviation (sigma) as a function of coupling strength and noise intensity
figure(1)
contourf(d_vec, del_vec, sig')
colorbar
xlabel('d')
ylabel('\delta_1')
zlabel('\sigma')
title('Mean Square Deviation')
set(gca, 'xscale', 'log', 'yscale', 'log')
clim([0.11 0.45]);

% Plot coefficient of variation (R) as a function of coupling strength and noise intensity
figure(2)
contourf(d_vec, del_vec, r')
colorbar
xlabel('d')
ylabel('\delta_1')
zlabel('R')
title('Coefficient of Variation')
set(gca, 'xscale', 'log', 'yscale', 'log')
clim([0.09 0.35]);
