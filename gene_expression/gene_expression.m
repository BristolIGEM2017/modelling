% Gene Expression Model
% Jeremie Joannes
% 29/08/2017
clear all; close all; clc;

%% Constants setup
% Nap transcription rates as polymerase rate x nucleotide count
k_DNA_Nap = 10; k_DNA_Nrf = 8;

% mRNA degradation rate
d_mRNA = 2;

% Translation rate -- IS IT SAME FOR BOTH NAP & NRF? --
k_mRNA = 8;

% Maturation rates
m_Nap = 20; m_Nrf = 20;

% Enzyme degradation rates
d_Nap = 10; d_Nrf = 10;

%% Simulation
% Set end time and timestep of simulation
t_end = 0.5; % s
dt = 0.0001; % s
tspan = 0:dt:t_end;

% Prepare initial conditions of
% DNA, mRNA_Nap, mRNA_Nrf, Nasc_Nap, Nasc_Nrf, Nap, Nrf]
init = [10, 0, 0, 0, 0, 0, 0];

% Run ODE solver, calling gene_expression_ODE.m
[t,C]=ode45(@(t,C) gene_expression_ODE(t,C,k_DNA_Nap,k_DNA_Nrf,d_mRNA,...
    k_mRNA,m_Nap,m_Nrf,d_Nap,d_Nrf),tspan,init);

%% Output
% Pass solver output into more readable variables
DNA = C(:,1);  % Not scrictly necessary
mRNA_Nap =  C(:,2);
mRNA_Nrf =  C(:,3);
Nasc_Nap =  C(:,4);
Nasc_Nrf =  C(:,5);
Nap =  C(:,6);
Nrf =  C(:,7);

% Plot outputs -- Nap & Nrf
figure
plot(t,Nap,'r',t,Nrf,'g','linewidth',2)
grid on
title('Gene Expression Model Results')
xlabel('Time [s]')
ylabel('Concentration [M]')
legend('Nap','Nrf')

