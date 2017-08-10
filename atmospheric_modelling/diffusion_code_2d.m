% Attempt at diffusion-only code in 2D
% Jeremie Joannes
% 09/08/2017
close all
clear all
clc

dx = 0.1;
nptsi = 101;
dy = 0.1;
nptsj = 101;

dt = 0.05;
nptst = 500;

a = .3; % diffusion coeff
amp = 3; % peak amplitude
u = .2; % x windspeed
v = .8; % y windspeed

% Initialise
f = zeros(nptsi,nptsj,nptst);
f(51,5,1) = amp; f(62,10,1) = 50*amp;

% Calcflate
for n = 1:nptst-1
    f(51,5,n) = amp; f(62,10,n) = amp;
    for i = 2:nptsi-1
        for j = 2:nptsj-1
            f(i,j,n+1) = f(i,j,n) + ...
                dt*a*(f(i+1,j,n)+f(i-1,j,n)-2*f(i,j,n))/dx + ...
                dt*a*(f(i,j+1,n)+f(i,j-1,n)-2*f(i,j,n))/dx;
            
            % If want advection:
            f(i,j,n+1) = f(i,j,n+1) - ...
                dt*(u*(f(i+1,j,n)-f(i-1,j,n))/(2*dx) + ...
                v*(f(i,j+1,n)-f(i,j-1,n))/(2*dy) );
        end
    end
end
        