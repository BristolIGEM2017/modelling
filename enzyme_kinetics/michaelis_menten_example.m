% Ecxample: Michaelis & Menten enzyme-substrate model

% close all
clear all
clc

k1 = 20;
km1 = 15;
k2 = 1;

Km = (k2+km1)/k1; % michaelis-menten constant

% Initial concentrations of S,E,C,P respectively
s_0 = 10; e_0 = 4; c_0 = 0; p_0 = 0;

Vmax = k2*e_0; %as predicted by QSSA

init = [s_0;e_0;c_0;p_0]; 
tspan = linspace(0,10,1000);

[t,C]=ode45(@(t,Y) exampleODE(t,Y,k1,km1,k2),tspan,init);

s = C(:,1); e = C(:,2); c = C(:,3); p = C(:,4);

V = Vmax.*s./(Km+s); % As predicted by QSSA

for i = 1:numel(t)
    dCdt(i,:) = exampleODE(t(i),[s(i),e(i),c(i),p(i)],k1,km1,k2);
end
dsdt = dCdt(:,1); dedt = dCdt(:,2); dcdt = dCdt(:,3); dpdt = dCdt(:,4);

figure
plot(t,s,t,e,t,c,t,p,'linewidth',2)
hold on
grid
title('Michaelis & Menten Enzyme-Substrate Complex Model')
xlabel('Time')
ylabel('Concentration')
legend('Substrate','Free enzyme','Complex','Product')

% Comparison with QSSA ---------------------------------
figure
plot(s,dpdt,s,V,'linewidth',2);
grid on
xlabel('[S]')
ylabel('V')
title('Comparison of exact ODE and QSSA - V/[S]')
legend('Actual dp/dt','M-M velocity')
axis([0 max(s) 0 max(V)])

[t_QSSA,C_QSSA]=ode45(@(t,C) example_QSSA_ODE(t,C,e_0,k1,km1,k2),tspan,[s_0;p_0]);
s_QSSA = C_QSSA(:,1); p_QSSA = C_QSSA(:,2);

figure
plot(t,p,t_QSSA,p_QSSA,'linewidth',2)
grid on
title('Comparison of exact ODE and QSSA - Time')
xlabel('Time')
ylabel('Product concentration')
legend('Exact ODE','QSSA')


