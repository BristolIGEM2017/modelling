% Attempt at diffusion-advection code in 2D
% Jeremie Joannes
% 09/08/2017
close all
clear all
clc

dx = 0.4; % m
dy = 0.4; % m

lenX = 30; % m
lenY = 30; % m

dt = 0.1; % s

a = .3; % diffusion coeff
u = 0.5; % x wind  speed in m/s
v = 0.2; % y windspeed in m/s

% Derived quantities
CFLx = 2*dt*a/dx^2;
CFLy = 2*dt*a/dy^2;
disp(['CFLx = ',num2str(CFLx),', CFLy = ',num2str(CFLy)]);
if CFLx > 1 || CFLy > 1
    error(['One/both diffusion CFL condition(s) not fulfilled - ',...
        'change dt, dx, dy and/or diffusivity coefficient']);
end

nptsi = lenX/dx + 4;
nptsj = lenY/dy + 4;

xpts = linspace(0,lenX,nptsi-4);
ypts = linspace(0,lenY,nptsj-4);

[Xmesh,Ymesh] = meshgrid(xpts,ypts);

% Read in sources and initialise f matrix
create_sources;

% Initialise time loop
E(1) = 1; n = 0; f_new = ones(nptsi,nptsj);

% Calculate
tic
while E > 1e-6
    % Update f
    n = n + 1;
    f = f_new;
    
    % Loop through i & j
    for i = 2:nptsi-1
        for j = 2:nptsj-1
            Diff = dt*a*(f(i+1,j)+f(i-1,j)-2*f(i,j))/dx^2 + ...
                dt*a*(f(i,j+1)+f(i,j-1)-2*f(i,j))/dy^2;
            
            Adv = dt*(u*(f(i+1,j)-f(i-1,j))/(2*dx) + ...
                v*(f(i,j+1)-f(i,j-1))/(2*dy));
            
            f_new(i,j) = f(i,j) + Diff - Adv;                
        end
    end
    
    % Impose source matrix for next time step 
    for srccounter = 1:numel(srcx)
        f_new(idx(srccounter),jdx(srccounter)) = srcval(srccounter);
    end
    
    % Update error
    E = max(max(abs(f-f_new))); Error_store(n) = E;
    disp(['Timestep = ', num2str(n),', Error = ', num2str(E)])
    
    % Instability condition break
    if E > 100
        error('Unstable scheme - review CFL or artificial diffusion');
        break
    end
end
toc

plot_diffusion_2d