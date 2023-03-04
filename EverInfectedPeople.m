% This script simulates the spread of an infectious disease using an SIR
% model with a modification for quarantined individuals.

% Set the simulation parameters
N = 2000;         % Total population size
beta = 0.4;         % Infection rate
gamma = 0.1;        % Recovery rate
%delta = 0.1;        % Quarantine rate
rho = 0.02;         % Death rate
%rho2 = 0.01;        % Death rate
I0 = 20;             % Initial number of infected individuals
S0 = N-I0;           % Initial number of susceptible individuals
R0 = 0;             % Initial number of recovered individuals
Q0 = 0;             % Initial number of quarantined individuals
D0 = 0;             % Initial number of death
C0 = 0;             % Initial cost
E0 = I0;             % Initial number of ever infected
T = 182;            % Simulation time in days
dt = 1/24;
clockmax = T/dt;

% Set up the arrays to store the results
S = zeros(1, T);
I = zeros(1, T);
R = zeros(1, T);
Q = zeros(1, T);
D = zeros(1, T);
P = zeros(1, T);
C = zeros(1, T);
E = zeros(1, T);
delta_arr = zeros(1, 1000);
d_arr = zeros(1,1000);
c_arr = zeros(1,1000);
e_arr = zeros(1,1000);
ptr = 1;

for delta = 0:0.001:1
    % Initialize the arrays with the initial values
    S(1) = S0;
    I(1) = I0;
    R(1) = R0;
    Q(1) = Q0;
    D(1) = D0;
    P(1) = S(1)+I(1)+R(1);
    C(1) = C0;
    E(1) = I0;
    % Run the simulation
    for t = 2:T
        % Calculate the number of new infections
        new_infections = 0;
        for i = 1:S(t-1)
            if rand < (beta*I(t-1)/P(t-1))
                new_infections = new_infections + 1;
            end
        end
    
        new_death_infections = 0;
        for i = 1:I(t-1)
            if rand < rho
                new_death_infections = new_death_infections+1;
            end
        end
    
        new_death_quarantines = 0;
        for i = 1:Q(t-1)
            if rand < rho
                new_death_quarantines = new_death_quarantines+1;
            end
        end

        % Calculate the number of individuals leaving the infected compartment
        new_recoveries = 0;
        for i = 1:I(t-1)
            if rand < gamma
                new_recoveries = new_recoveries + 1;
            end
        end
    
        % Calculate the number of individuals leaving the quarantined compartment
        new_dequarantines = 0;
        for i = 1:Q(t-1)
            if rand < gamma
                new_dequarantines = new_dequarantines + 1;
            end
        end
    
        % Calculate the number of new quarantines
        new_quarantines = 0;
        for i = 1:I(t-1)
            if rand < delta
                new_quarantines = new_quarantines + 1;
            end
        end

        % Update the compartment sizes
        S(t) = S(t-1) - new_infections;
        I(t) = I(t-1) + new_infections - new_recoveries - new_death_infections - new_quarantines;
        R(t) = R(t-1) + new_recoveries + new_dequarantines;
        Q(t) = Q(t-1) + new_quarantines - new_dequarantines - new_death_quarantines;
        P(t) = S(t-1) + I(t-1) + R(t-1);
        D(t) = D(t-1) + new_death_quarantines + new_death_infections;
        C(t) = C(t-1) + Q(t-1);
        E(t) = E(t-1) + I(t-1);
   

    end
    delta_arr(ptr) = delta;
    d_arr(ptr) = D(T);
    c_arr(ptr) = C(T);
    e_arr(ptr) = E(T);    
    %delta_arr(ptr)
    %d_arr(ptr)
    ptr = ptr+1;
end

figure;

plot(delta_arr, e_arr, 'b-', 'LineWidth', 2);
xlabel('Quarantine rate');
ylabel('Total Infected Days');
title('Benefit and Cost of Quarantine rate');
hold on;
plot(delta_arr, c_arr, 'r-', 'LineWidth', 2);
xlabel('Quarantine rate');
ylabel('Days');
legend('Total Infected Days', 'Total Quarantine Days');

   % %ylabel('Quarantine People');
%title('Cost of Quarantine rate');


%{
% Plot the result
subplot(2,1,1),...
    plot(delta_arr, e_arr, 'b-', 'LineWidth', 2);
    xlabel('Quarantine rate');
    ylabel('Total Infected Days');
    title('Benefit of Quarantine rate');

subplot(2,1,2),...
    plot(delta_arr, c_arr, 'r-', 'LineWidth', 2);
    xlabel('Quarantine rate');
    ylabel('Total Quarantine Days');
   %%ylabel('Quarantine People');
    title('Cost of Quarantine rate');
%}

