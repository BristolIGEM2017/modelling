%% Pollution Display
figure(1)
fontsz = 11;
to_plot = f(3:nptsi-2,3:nptsj-2);
contour(Xmesh',Ymesh',to_plot,20,'fill','on');
caxis([0 max(max(f(:,:,1)))])
set(gca,'fontsize',fontsz);

colormap(flipud(gray))

c = colorbar;
set(c,'fontsize',fontsz);
ylabel(c,'Concentration [\mug/m^3]')
ylabel('y-coordinate [m]')
xlabel('x-coordinate [m]')

titlestr = ['Pollution map at timestep ',num2str(n),...
    ' with error ', num2str(E)];
title(titlestr)

grid on
axis equal

%% Convergence History
figure(2)
semilogy(1:n,Error_store,'k-')
set(gca,'fontsize',fontsz);
title('Convergence History')
ylabel('Max grid error')
xlabel('Timestep')
grid on 