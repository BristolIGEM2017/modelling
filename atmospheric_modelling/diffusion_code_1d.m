% Attempt at diffusion-only code in 1D
% Jeremie Joannes
% 09/08/2017

dx = 0.1;
nptsi = 101;

dt = 0.05;
nptst = 500;

a = .2; % diffusion coeff
v = 0.4; % transport velocity

% Initialise
u = zeros(nptsi,nptst);
u(20,1) = 1; u(30,1) = 2;

% Calculate
for n = 1:nptst-1
    for i = 2:nptsi-1
        
        u(i,n+1) = u(i,n) + dt*a*(u(i+1,n)+u(i-1,n)-2*u(i,n))/dx;
        % Add convection component:
        u(i,n+1) = u(i,n+1) - dt*v*(u(i+1,n)-u(i-1,n))/(2*dx);
    end
end

for n = 1:nptst
    figure(1)
    plot(2:nptsi-1,u(2:nptsi-1,n),'.-')
    axis([1 101 -0.2 2])
    pause(0.05)
end
        