% Attempt at diffusion-advection code in 2D
% Jeremie Joannes
% 09/08/2017
close all
clear all
clc

dx = 0.5; % m
dy = 0.5; % m

lenX = 30; % m
lenY = 30; % m

dt = 0.05; % s
nptst = 1000;

a = 0.3; % diffusion coeff
u = -.02; % x wind  speed in m/s
v = -.5; % y windspeed in m/s

% Derived quantities
nptsi = lenX/dx + 4;
nptsj = lenY/dy + 4;

xpts = linspace(0,lenX,nptsi-4);
ypts = linspace(0,lenY,nptsj-4);

[Xmesh,Ymesh] = meshgrid(xpts,ypts);

% Read in sources and initialise f matrix
create_sources;

% Calculate
tic
for n = 1:nptst-1
    for i = 2:nptsi-1
        for j = 2:nptsj-1
            f(i,j,n+1) = f(i,j,n) + ...
                dt*a*(f(i+1,j,n)+f(i-1,j,n)-2*f(i,j,n))/dx^2 + ...
                dt*a*(f(i,j+1,n)+f(i,j-1,n)-2*f(i,j,n))/dy^2;
            
            % If also want advection:
            f(i,j,n+1) = f(i,j,n+1) - ...
                dt*(u*(f(i+1,j,n)-f(i-1,j,n))/(2*dx) + ...
                v*(f(i,j+1,n)-f(i,j-1,n))/(2*dy));
        end
    end
    for srccounter = 1:numel(srcx)
        % Impose source matrix for next time step
        f(idx(srccounter),jdx(srccounter),n+1) = srcval(srccounter);
    end
end
toc

plot_diffusion_2d