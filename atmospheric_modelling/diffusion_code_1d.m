% Attempt at diffusion-only code in 1D
% Jeremie Joannes
% 09/08/2017
clear all
close all
clc

dx = 1;
nptsi = 101;

dt = 1;

a = 0.15; % diffusion coeff
v = 0.3; % transport velocity

% CFL for the diffusion equation
CFL = 2*dt*a/dx^2;
if CFL > 1
    error(['CFL of ',num2str(CFL),...
        ' too large, reduce dt, diffusivity or increase dx']);
else
    disp(['CFL of the diffusion scheme: ',num2str(CFL)])
end

% Initialise
u = zeros(nptsi); u_new = ones(nptsi);
u(20) = 1; u(30) = 2;

% Calculate
E = 1; n = 0;
while E > 1e-7 % Tolerance
    % Update everything
    u = u_new;
    n = n + 1;

    % Start the x loop
    for i = 2:nptsi-1
        diff = dt*a*(u(i+1)+u(i-1)-2*u(i))/dx^2;
        adv = dt*v*(u(i+1)-u(i-1))/(2*dx);
        u_new(i) = u(i) + diff - adv;
    end
    % Re-impose conditions
    u_new(20) = 1; u_new(30) = 2;
    
    % Calculate & display current max error
    E = max(abs(u(:)-u_new(:))); Error_store(n) = E;
    disp(['Timestep = ', num2str(n),', Error = ', num2str(E)])
    
    % Break condition in case unstable
    if E > 10
        error('Unstable scheme - review CFL or artificial diffusion');
        break
    end
end

%% Plots
figure(1)
semilogy(1:n,Error_store,'k-')
title('Convergence History')
ylabel('Max grid error')
xlabel('Timestep')
grid on

figure(2)
plot(2:nptsi-1,u(2:nptsi-1),'.-')
axis([1 101 -0.2 2])
titlestr = ['Max grid error at end timestep ', num2str(n), ' = ',num2str(E)];
title(titlestr)  