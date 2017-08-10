for n = 1:nptst
    figure(1)
    to_plot = f(:,:,n);
    contour(to_plot,20,'fill','on')
    colormap(flipud(gray))
    caxis([0 max(max(f(:,:,1)))])
    axis equal
    pause(0.05)
end