%clear all; close all
% The following code computes the mean phase coherence for the L-O system

% time
dt = 0.01;       % timestep
T = 200;         % duration of simulation
t0 = 50/dt;      % transient time (index of transient)
t = (0:dt:T);    % time domain

% coupling
d = 0.1;       % common coupling strength
d21 = d;       % coupling between Osc. 1 & 2
d31 = d;       % coupling between Osc. 1 & 3
d32 = d;       % coupling between Osc. 2 & 3

num_trials = 25;               % number of trials
del_vec = logspace(-3,0.6,20); % vector of noise intensities

for j = 1:length(del_vec) % iterate over noise intensities

    % noise intensity
    del_1 = del_vec(j);  % driving noise intensity (Osc.1)
    del_2 = 0.01;        % noise intensity Osc.2
    del_3 = 0.01;        % noise intensity Osc.3

    for i = 1:num_trials % iterate over trials

        % SIMULATE MODEL
        [X] = L_O_Motif_simulation(del_1,del_2,del_3,d21,d31,d32,dt,T);

        % SMOOTH TIMESERIES (after removing transient time)
        x1 = smoothdata(X(1,t0:end),'gaussian',100);
        y1 = smoothdata(X(2,t0:end),'gaussian',100);
        x2 = smoothdata(X(3,t0:end),'gaussian',100);
        y2 = smoothdata(X(4,t0:end),'gaussian',100);
        x3 = smoothdata(X(5,t0:end),'gaussian',100);
        y3 = smoothdata(X(6,t0:end),'gaussian',100);

        % CALCULATE PHASE OF OSCILLATORS
        p1 = atan2(y1,x1);
        p2 = atan2(y2,x2);
        p3 = atan2(y3,x3);

        % CALCULATE PHASE DIFFERENCE OF OSCILLATORS
        %pd12 = 0.5*(p1-p2); % phase difference osc. 1 & 2
        pd13 = 0.5*(p1-p3); % phase difference osc. 1 & 3
        %pd23 = 0.5*(p2-p3); % phase difference osc. 2 & 3

        % AVERAGE ABSOLUTE PHASE DIFFERENCE OVER TIME (in absolute value)
        %avg_pd12_temp(i,j) = mean(abs(pd12));   % ' ' osc. 1 & 2
        avg_pd13_temp(i,j) = mean(abs(pd13));   % ' ' osc. 1 & 3
        %avg_pd23_temp(i,j) = mean(abs(pd23));   % ' ' osc. 2 & 3

        % CALCULATE PHASE COHERENCE
        R_temp(i,j) = sqrt(mean(cos(pd13))^2 + mean(sin(pd13))^2); % phase coherence

    end % end iteration over trials

    % AVERAGE PHASE DIFFERENCE OVER TRIALS (in absolute value)
    %avg_pd12 = mean(avg_pd12_temp);   % ' ' osc. 1 & 2
    avg_pd13 = mean(avg_pd13_temp);    % ' ' osc. 1 & 3
    %avg_pd23 = mean(avg_pd23_temp);   % ' ' osc. 2 & 3
    R = mean(R_temp); % mean phase coherence

end % end of iteration over noise

% PLOT MEAN ABSOLUTE PHASE DIFFERENCE
figure(1)
semilogx(del_vec,avg_pd13)
ylabel('|\phi_1-\phi_3|')
xlabel('\delta_1')
title('Mean Absolute Phase Difference')

% PLOT MEAN PHASE COHERENCE
figure(2)
semilogx(del_vec,R)
ylabel('\gamma')
xlabel('\delta_1')
title('Mean Phase Coherence')
