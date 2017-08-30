% Mapping code to convert road shapefiles into cartesian coordinates
% Albert Wigmore
% 30/08/17

close all
clear all

% Names of different road types
roadTypes = {'A Road', 'A Road, Collapsed Dual Carriageway', ...
             'B Road', 'B Road, Collapsed Dual Carriageway', ...
             'Local Street', ...
             'Minor Road', 'Minor Road, Collapsed Dual Carriageway', ...
             'Motorway', 'Motorway, Collapsed Dual Carriageway', ...
             'Pedestrianised Street', ...
             'Primary Road', 'Primary Road, Collapsed Dual Carriageway', ...
             'Private Road Publicly Accessible'};
         
% Source values corresponding to road types
roadVals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
mapObj = containers.Map(roadTypes, roadVals);

% Read in 'Road' shapefile
S = shaperead('mapping/Road.shp');

% Loop through all the roads and create sources
output = [];
h = waitbar(0, 'Please wait, creating sources...');
for i = 1:size(S, 1)
    for j = 1:size(S(i).X, 2)-1
        [x, y] = bresenham(S(i).X(j), S(i).X(j + 1), ...
                           S(i).Y(j), S(i).Y(j + 1));
        if not(isnan(x))
            output = [output; x, y, ones(size(x)) * mapObj(S(i).classifica)];
        end
    end
    waitbar(i / size(S, 1))
end
close(h);

% Write data to CSV file
csvwrite('sources.dat', output)
