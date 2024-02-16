function [X] = L_O_Motif_simulation(del_1, del_2, del_3, d21, d31, d32, dt, T)
% L_O_Motif2 simulates a 3-node motif (version 2) using a Lorenz-Oscillator (L-O) system
%   del_i is the noise intensity corresponding to oscillator i
%   dji is the coupling strength between oscillator i & j
%   dt is the time step
%   T is the duration of simulation (i.e. (0:dt:T))

% TIME PARAMETERS
t = (0:dt:T);   % Time domain

% MODEL PARAMETERS
lam_0 = -0.1;   % Lambda
alp = -0.2;     % Alpha
gam = -0.2;     % Gamma
om0 = 2;        % Omega_0
om1 = 0;        % Omega_1

% DEPENDENT VARIABLES/MODEL OUTPUT (X = [x1, y1, x2, y2, x3, y3])
X = zeros(6, length(t));

% E-M METHOD/MODEL SIMULATION
% Set initial conditions
X(:, 1) = [0.0090, 0.0206, -0.0923, 0.1282, 0.0089, 0.0064]; % Small (arbitrary) initial conditions

% Generate noise vector for each oscillator
noise = randn(3, length(t));

% Iterate over time steps
for k = 2:length(t)
    % Calculate squared radii for each oscillator
    rsq = [X(1, k-1)^2 + X(2, k-1)^2;
        X(3, k-1)^2 + X(4, k-1)^2;
        X(5, k-1)^2 + X(6, k-1)^2]; % (r_i)^2 = (x_i)^2 + (y_i)^2

    % Update state variables for each oscillator
    X(:, k) = [
        % Oscillator 1
        X(1, k-1) + ((lam_0 + alp*rsq(1) + gam*rsq(1)^2)*X(1, k-1) - ...
        (om0 + om1*rsq(1))*X(2, k-1) + 0*(X(3, k-1)- X(1, k-1)))*dt  + ...
        del_1*sqrt(dt)*noise(1, k-1); % x1

        X(2, k-1) + ((om0 + om1*rsq(1))*X(1, k-1) + (lam_0 + alp*rsq(1) + ...
        gam*rsq(1)^2)*X(2, k-1) + 0*(X(4, k-1)- X(2, k-1)))*dt;   % y1

        % Oscillator 2
        X(3, k-1) + ((lam_0 + alp*rsq(2) + gam*rsq(2)^2)*X(3, k-1) - ...
        (om0 + om1*rsq(2))*X(4, k-1) + d21*(X(1, k-1)- X(3, k-1)))*dt + ...
        del_2*sqrt(dt)*noise(2, k-1);  % x2

        X(4, k-1) + ((om0 + om1*rsq(2))*X(3, k-1) + (lam_0 + alp*rsq(2) ...
        + gam*rsq(2)^2)*X(4, k-1) + d21*(X(2, k-1)- X(4, k-1)))*dt; % y2

        % Oscillator 3
        X(5, k-1) + ((lam_0 + alp*rsq(3) + gam*rsq(3)^2)*X(5, k-1) - ...
        (om0 + om1*rsq(3))*X(6, k-1) + d31*(X(1, k-1)- X(5, k-1)) -...
        d32*(X(3, k-1)- X(5, k-1)))*dt + del_3*sqrt(dt)*noise(3, k-1); % x3

        X(6, k-1) + ((om0 + om1*rsq(3))*X(5, k-1) + (lam_0 + alp*rsq(3) + ...
        gam*rsq(3)^2)*X(6, k-1) + d31*(X(2, k-1)- X(6, k-1)) ...
        - d32*(X(4, k-1)- X(6, k-1)))*dt; % y3
        ];
end

end % End of function
