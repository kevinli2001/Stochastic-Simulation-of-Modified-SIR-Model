clear all;
% This script simulates the spread of an infectious disease using an SIR
% model with a modification for quarantined individuals.

% Set the simulation parameters
N = 2000;         % Total population size
beta = 0.4;         % Infection rate
gamma = 0.1;        % Recovery rate
delta = 0.1;        % Quarantine rate
rho = 0.02;         % Death rate
%rho2 = 0.01;        % Death rate
I0 = 20;             % Initial number of infected individuals
S0 = N-I0;           % Initial number of susceptible individuals
R0 = 0;             % Initial number of recovered individuals
Q0 = 0;             % Initial number of quarantined individuals
D0 = 0;             % Initial number of death
C0 = 0;             % Initial cost
T = 182;            % Simulation time in days
dt = 1/24;
clockmax = T/dt;


% Set up the arrays to store the results
S = zeros(1, clockmax);
I = zeros(1, clockmax);
R = zeros(1, clockmax);
Q = zeros(1, clockmax);
D = zeros(1, clockmax);
P = zeros(1, clockmax);
C = zeros(1, clockmax);

% Initialize the arrays with the initial values
S(1) = S0;
I(1) = I0;
R(1) = R0;
Q(1) = Q0;
D(1) = D0;
P(1) = S(1)+I(1)+R(1);
C(1) = C0;

% Run the simulation
for t = 2:clockmax
    % Calculate the number of new infections
    new_infections = beta*S(t-1)*I(t-1)/P(t-1)*dt;

    % Calculate the number of death
    new_death = (rho*I(t-1) + rho*Q(t-1))*dt;
    
    % Calculate the number of individuals leaving the infected compartment
    new_recoveries = gamma*I(t-1)*dt; %-rho*I(t-1)-rho2*I(t-1);
    
    % Calculate the number of individuals leaving the quarantined compartment
    new_quarantines = delta*I(t-1)*dt;
    new_dequarantines = gamma*Q(t-1)*dt; % number of quarantined individuals who recover
    
    % Update the compartment sizes
    S(t) = S(t-1) - new_infections;
    I(t) = I(t-1) + new_infections - new_recoveries - new_quarantines - rho*I(t-1)*dt;
    R(t) = R(t-1) + new_recoveries + new_dequarantines;
    Q(t) = Q(t-1) + new_quarantines - new_dequarantines - rho*Q(t-1)*dt;
    D(t) = D(t-1) + new_death;
    P(t) = S(t-1) + I(t-1) + R(t-1);
    C(t) = C(t-1) + Q(t-1) - rho*Q(t-1);
end

%N = S + I + R;

% Plot the results
figure;
plot(S, 'b', 'LineWidth', 2);
hold on;
plot(I, 'r', 'LineWidth', 2);
plot(R, 'g', 'LineWidth', 2);
plot(Q, 'm', 'LineWidth', 2);
plot(D, 'c', 'LineWidth', 2);
plot(P, 'k', 'LineWidth', 2);
%plot(C, 'y', 'LineWidth', 2);
%legend('Susceptible', 'Infected', 'Recovered');
legend('Susceptible', 'Infected', 'Recovered', 'Quarantined', 'Death', 'Population');
xlabel('Time (Hours)');
ylabel('Number of individuals');
title('Spread of an Infectious Diseas');
%title('Spread of an Infectious Disease with Quarantine');
