function [] = enzyme_kinetics(Kcat_Nap,Kd_Nap,Km_Nap,Kcat_Nrf,Kd_Nrf,Km_Nrf)
    % EITHER: define rate constants (in terms of grouped constants)
    k1 = Kcat_Nap/(Km_Nap-Kd_Nap);
    k_1 = Kd_Nap*k1;
    k2 = Kcat_Nap;

    k3 = Kcat_Nrf/(Km_Nrf-Kd_Nrf);
    k_3 = Kd_Nrf*k3;
    k4 = Kcat_Nrf;

    % OR: define rate constants (in terms of individual paths)
    % k1 = 5;
    % k_1 = 3;
    % k2 = 2;
    % k3 = 5;
    % k_3 = 3;
    % k4 = 2;

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
    init = [NO3_0; Nap_0; NapNO3_0; NO2_0; Nrf_0; NrfNO2_0; NH4_0];
    tspan = 0:dt:t_end;

    % Run ODE solver, calling enzyme_ODE.m
    [t,C]=ode45(@(t,C) enzyme_ODE(t, C, k1, k_1, k2, k3, k_3, k4),tspan,init);

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
    figure(1)
    plot(t,NO3,'r',t,NO2,'y',t,NH4,'g','linewidth',0.5)
    hold on
    grid on
    title('Michaelis & Menten model showing nitrate, nitrite, ammonia')
    xlabel('Time [s]')
    ylabel('Concentration [M]')
    legend('NO_3^-','NO_2^-','NH_4^+')

    % Plot outputs -- enzymes & complexes
    figure(2)
    plot(t,Nap,'b',t,NapNO3,'r',t,Nrf,'g',t,NrfNO2,'y','linewidth',0.5)
    hold on
    grid on
    title('Michaelis & Menten model showing enzymes and complexes')
    xlabel('Time [s]')
    ylabel('Concentration [M]')
    legend('Free Nap','Nap \cdot NO_3^-','Free Nrf','Nrf \cdot NO_2^-')

    % Plot outputs -- reaction velocity
    figure(3)
    plot(t(1:end-1),V,'b','linewidth',0.5)
    hold on
    grid on
    title('Reaction velocity')
    xlabel('Time [s]')
    ylabel('Velocity [M/s]')
end