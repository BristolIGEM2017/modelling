fontsz = 16;
no_of_disp_times = 20; % Choose how many timepoints you want to display

for tpt = 1:no_of_disp_times
    figure(1)

    to_plot = f(3:nptsi-2,3:nptsj-2,tpt*nptst/no_of_disp_times);
    contour(Xmesh',Ymesh',to_plot,20,'fill','on');
    caxis([0 max(max(f(:,:,1)))])
    set(gca,'fontsize',fontsz);
    
    colormap(flipud(gray))
    
    c = colorbar;
    set(c,'fontsize',fontsz);
    ylabel(c,'Concentration [\mug/m^3]')
    ylabel('y-coordinate [m]')
    xlabel('x-coordinate [m]')
        
    titlestr = ['Pollution map at time t = ', num2str(tpt*nptst/no_of_disp_times*dt), ' s'];
    title(titlestr)
    
    grid on
    axis equal
    
    pause(.01);
end
