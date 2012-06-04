function [ ] = PlotMultiPanelSample( Var, Lon, Lat, Ccn, Sst, OutFile )
%PlotMultiPanelSample Create a multipanel plot of the sample data
%   This function will create a multi-panel plot of the data given in Var.
%   Var is a 3D array (n,x,y) where n equals the length of Ccn and Sst. Lon
%   holds the longitude values (x) and Lat holds the latitude values (y).

% For now assume that n is 4 --> 2x2 subplots.

Fig = figure;

% Change the zero value to white
Cmap = colormap;
Cmap(1,:) = [ 1 1 1 ];
colormap(Cmap);

for i = 1:4
    subplot(2,2,i);
    contourf(Lon,Lat,squeeze(Var(i,:,:)));
    Xlab = sprintf('CCN: %d/cc, SST: %d',Ccn(i),Sst(i));
    xlabel(Xlab);
    colorbar;
end

saveas(Fig,OutFile);

colormap('default');
close(Fig);


end

