function [colourplot, scanplot] = showscan(datalist,startpoint,endpoint,np)

% Construct a scan from 2D data and do two plots:
%  the scan and the original map with the graphical representation of the scan path
%
% P. Steffens, 03/2008


startpoint = startpoint(:)'; 
endpoint   = endpoint(:)';

% construct the scan points
scanpath = repmat(startpoint,np,1) + (0:(np-1))'/np  * (endpoint-startpoint);

% Ensure qx,qy-coordinates for the data points
datalist = coordtransform(datalist,'qplane');

% Obtain values and errors by interpolation
[z,dz]  = linearinterpolation (datalist, [scanpath(:,1), scanpath(:,2)]);

% Do the color plot
colourplot = fcplot(datalist);
hold on
plot(scanpath(:,1), scanpath(:,2), '-ok');
hold off

% Do the scan plot
scanplot = figure ;
[m,i] = max( max(scanpath,[],1) ); %i is No. of column with largest variation; take this as x-axis
errorbar(scanpath(:,i),z,dz,'ob');
xlabels= {'Qx','Qy'};
xlabel([xlabels{i} '(' char(197) '^{-1})']);