clear all
close all
clc

Kcat_Nap = linspace(1.9,2.1,3);
Kd_Nap = linspace(0.5,0.7,3);
Km_Nap = linspace(0.8,1.2,3);

Kcat_Nrf = linspace(1.9,2.1,3);
Kd_Nrf = linspace(0.5,0.7,3);
Km_Nrf = linspace(0.8,1.2,3);

for Kcat_Nap = Kcat_Nap(1:end)
    enzyme_kinetics(Kcat_Nap,mean(Kd_Nap),mean(Km_Nap)...
        ,mean(Kcat_Nrf),mean(Kd_Nrf),mean(Km_Nrf));
end
for Kd_Nap = Kd_Nap(1:end)
    enzyme_kinetics(mean(Kcat_Nap),(Kd_Nap),mean(Km_Nap)...
        ,mean(Kcat_Nrf),mean(Kd_Nrf),mean(Km_Nrf));
end
for Km_Nap = Km_Nap(1:end)
    enzyme_kinetics(mean(Kcat_Nap),mean(Kd_Nap),(Km_Nap)...
        ,mean(Kcat_Nrf),mean(Kd_Nrf),mean(Km_Nrf));
end

for Kcat_Nrf = Kcat_Nrf(1:end)
    enzyme_kinetics(mean(Kcat_Nap),mean(Kd_Nap),mean(Km_Nap)...
        ,(Kcat_Nrf),mean(Kd_Nrf),mean(Km_Nrf));
end
for Kd_Nrf = Kd_Nrf(1:end)
    enzyme_kinetics(mean(Kcat_Nap),mean(Kd_Nap),mean(Km_Nap)...
        ,mean(Kcat_Nrf),(Kd_Nrf),mean(Km_Nrf));
end
for Km_Nrf = Km_Nrf(1:end)
    enzyme_kinetics(mean(Kcat_Nap),mean(Kd_Nap),(Km_Nap)...
        ,mean(Kcat_Nrf),mean(Kd_Nrf),(Km_Nrf));
end