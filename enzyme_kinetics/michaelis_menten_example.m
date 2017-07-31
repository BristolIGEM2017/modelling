% Ecxample: Michaelis & Menten enzyme-substrate model

close all
clear all
clc

k1 = 6;
km1 = 1;
k2 = 10;

init = [10;4;0;0]; %Concentrations of S,E,C,P respectively
tspan = linspace(0,1,1000);

[t,C]=ode45(@(t,Y) exampleODE(t,Y,k1,km1,k2),tspan,init);

s = C(:,1); e = C(:,2); c = C(:,3); p = C(:,4);
figure
plot(t,s,t,e,t,c,t,p,'linewidth',2)
hold
grid
title('Michaelis & Menten Enzyme-Substrate Complex Model')
xlabel('Time')
ylabel('Concentration')
legend('Substrate','Free enzyme','Complex','Product')
hold off
