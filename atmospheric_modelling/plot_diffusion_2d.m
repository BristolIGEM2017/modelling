

% for n = 1:nptst
%     figure(1)
%     to_plot = f(3:nptsi-2,3:nptsj-2,n);
%     contour(xpts,ypts,to_plot,20,'fill','on')
%     colormap(flipud(gray))
%     caxis([0 max(max(f(:,:,1)))])
%     axis equal
%     pause(0.05)
% end

    figure(1)
    to_plot = f(3:nptsi-2,3:nptsj-2,end);
    contour(Xmesh',Ymesh',to_plot,20,'fill','on')
    colormap(flipud(gray))
    caxis([0 max(max(f(:,:,1)))])
    axis equal

