% Load source data
M = csvread('sources.csv',1);
srcx = M(:,1); srcy = M(:,2); srcint = M(:,3);

S = zeros(nptsi,nptsj,nptst); % Initialise source vector

for counter = 1:numel(srcx)
    tmp = abs(xpts-srcx(counter));
    [idx(counter),idx(counter)] = min(tmp); %index of closest value
    idx(counter) = idx(counter) + 2; %closest value in terms of indices incl. 2 halo cells

    tmp = abs(ypts-srcy(counter));
    [jdx(counter) ,jdx(counter)] = min(tmp); %index of closest value
    jdx(counter) = jdx(counter) + 2; %closest value in terms of indices incl. 2 halo cells

    S(idx(counter),jdx(counter),1) = srcint(counter);% Create source matrix
end 
