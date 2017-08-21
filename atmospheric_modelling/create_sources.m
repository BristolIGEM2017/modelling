% Load source data
M = csvread('sources.csv',1);
srcx = M(:,1); srcy = M(:,2); srcval = M(:,3);

f = zeros(nptsi,nptsj); % Initialise source vector

for srccounter = 1:numel(srcx)
    tmp = abs(xpts-srcx(srccounter));
    [idx(srccounter),idx(srccounter)] = min(tmp); %index of closest value
    idx(srccounter) = idx(srccounter) + 2; %closest value in terms of indices incl. 2 halo cells

    tmp = abs(ypts-srcy(srccounter));
    [jdx(srccounter) ,jdx(srccounter)] = min(tmp); %index of closest value
    jdx(srccounter) = jdx(srccounter) + 2; %closest value in terms of indices incl. 2 halo cells

    f(idx(srccounter),jdx(srccounter)) = srcval(srccounter);% Create source matrix
end 
