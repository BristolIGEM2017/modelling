function dCdt = enzyme_ODE(t,C,k1,km1,k2,k3,km3,k4)
    % Define each concentration for ease of use
    NO3 = C(1);
    Nap =  C(2);
    NapNO3 =  C(3);
    NO2 =  C(4);
    Nrf =  C(5);
    NrfNO2 =  C(6);
    NH4 =  C(7); % Not scrictly necessary

    % Define each species' derivative
    dNO3_dt = km1*NapNO3 - k1*Nap*NO3;
    dNap_dt = km1*NapNO3 - k1*Nap*NO3 + k2*NapNO3;
    dNapNO3_dt = -km1*NapNO3 + k1*Nap*NO3 - + k2*NapNO3;
    dNO2_dt = k2*NapNO3 + km3*NrfNO2  - k3*Nrf*NO2;
    dNrf_dt = km3*NrfNO2  - k3*Nrf*NO2 + k4*NrfNO2;
    dNrfNO2_dt = k3*Nrf*NO2 - km3*NrfNO2 - k4*NrfNO2;
    dNH4_dt = k4*NrfNO2;

    % Lump together into one big d/dt vector
    dCdt = [dNO3_dt; dNap_dt; dNapNO3_dt; dNO2_dt; dNrf_dt; dNrfNO2_dt; dNH4_dt];
end