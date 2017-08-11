% Attempt at diffusion-advection code in 2D
% Jeremie Joannes
% 09/08/2017
close all
clear all
clc

dx = 0.1; % m
dy = 0.1; % m

lenX = 25; % m
lenY = 25; % m

dt = 0.05; % s
nptst = 1000;

a = 0.3; % diffusion coeff
src = 100; % source amplitude in ppm
u = .2; % x windspeed in m/s
v = .9; % y windspeed in m/s

% Derived quantities
nptsi = lenX/dx + 4;
nptsj = lenY/dy + 4;

xpts = linspace(0,lenX,nptsi-4);
ypts = linspace(0,lenY,nptsj-4);

[Xmesh,Ymesh] = meshgrid(xpts,ypts);

% Read in sources
create_sources;

% Initialise
f = S; % Equal to the source matrix

% Calculate
tic
for n = 1:nptst-1
    for counter = 1:numel(srcx)
        f(idx(counter),jdx(counter),n) = srcint(counter); % Impose source matrix
    end
    for i = 2:nptsi-1
        for j = 2:nptsj-1
            f(i,j,n+1) = f(i,j,n) + ...
                dt*a*(f(i+1,j,n)+f(i-1,j,n)-2*f(i,j,n))/dx + ...
                dt*a*(f(i,j+1,n)+f(i,j-1,n)-2*f(i,j,n))/dx;
            
            % If want advection:
            f(i,j,n+1) = f(i,j,n+1) - ...
                dt*(u*(f(i+1,j,n)-f(i-1,j,n))/(2*dx) + ...
                v*(f(i,j+1,n)-f(i,j-1,n))/(2*dy));
        end
    end
end
toc

plot_diffusion_2d