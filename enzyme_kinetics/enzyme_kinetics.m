% 
close all
clear all
clc

% Define rate constants
k1 = 5;
km1 = 3;
k2 = 2;
k3 = 5;
km3 = 3;
k4 = 2;

% Initial concentrations of each species
NO3_0 = 10;
Nap_0 =  3;
NapNO3_0 =  0;
NO2_0 =  5;
Nrf_0 =  3;
NrfNO2_0 =  0;
NH4_0 =  0;

% Set end time and timestep of simulation
t_end = 5; % s
dt = 0.001; % s

% Prepare initial conditions
init = [NO3_0;Nap_0;NapNO3_0;NO2_0;Nrf_0;NrfNO2_0;NH4_0];
tspan = 0:dt:t_end;

% Run ODE solver, calling enzyme_ODE.m
[t,C]=ode45(@(t,C) enzyme_ODE(t,C,k1,km1,k2,k3,km3,k4),tspan,init);

% Pass solver output into more readable variables
NO3 = C(:,1);
Nap =  C(:,2);
NapNO3 =  C(:,3);
NO2 =  C(:,4);
Nrf =  C(:,5);
NrfNO2 =  C(:,6);
NH4 =  C(:,7); % Not scrictly necessary

% Calculate reaction velocity
V = diff(NH4)/dt;

% Plot outputs -- reactants & products
figure
plot(t,NO3,t,NO2,t,NH4,'linewidth',2)
hold on
grid
title('Michaelis & Menten model showing nitrate, nitrite, ammonia')
xlabel('Time [s]')
ylabel('Concentration [M]')
legend('NO_3^-','NO_2^-','NH_4^+')

% Plot outputs -- enzymes & complexes
figure
plot(t,Nap,t,NapNO3,t,Nrf,t,NrfNO2,'linewidth',2)
hold on
grid
title('Michaelis & Menten model showing enzymes and complexes')
xlabel('Time [s]')
ylabel('Concentration [M]')
legend('Free Nap','Nap \cdot NO_3^-','Free Nrf','Nrf \cdot NO_2^-')

% Plot outputs -- reaction velocity
figure
plot(t(1:end-1),V,'b','linewidth',2)
grid
title('Reaction velocity')
xlabel('Time [s]')
ylabel('Velocity [M/s]')