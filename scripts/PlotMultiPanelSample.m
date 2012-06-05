function [ ] = PlotMultiPanelSample( Var, Lon, Lat, Ccn, Sst, Cmap, Clevs, Ptitle, OutFile )
%PlotMultiPanelSample Create a multipanel plot of the sample data
%   This function will create a multi-panel plot of the data given in Var.
%   Var is a 3D array (n,x,y) where n equals the length of Ccn and Sst. Lon
%   holds the longitude values (x) and Lat holds the latitude values (y).

% For now assume that n is 4 --> 2x2 subplots.

Fig = figure;

% Set the colormap
colormap(Cmap);

% Set up portion of subplot regions for 4 panels
% with a title above and a colorbar below the 4 panels.
% The idea is to make the field plots cover more subplot
% spaces than the title and colorbar so that field plots can be
% larger. Define a 6x2 subplot scheme, use 2 vertical subplot
% regions for the panels and 1 veritcal subplot region for each
% of the title and colorbar. Subplot arranges looks like:
%
%   1  2
%   3  4
%   5  6
%   7  8
%   9 10
%  11 12
%
% where the first plot goes in the 2x1 region covered by subplots
% indices: 3 and 5. Second plot goes in the 2x1 region with
% indices:  4 and 6. Etc.
%
% The title goes in the 1x2 region at the top (indices 1 and 2)
%
% The colorbar goes in the 1x2 region at the bottom (indices
% 11 and 12).
%

Sregion = [ 3 5;
4 6;
7 9;
8 10 ];

Fsize = 18;

% Place the text in the center of the top row
subplot(6,2, [ 1 2 ]);
axis off;
text(0.5,0.5,Ptitle,'HorizontalAlignment', 'center', ...
  'VerticalAlignment', 'middle', 'FontSize', 24);

for i = 1:4
    subplot(6,2,squeeze(Sregion(i,:)));
    contourf(Lon,Lat,squeeze(Var(i,:,:)),Clevs), shading flat;
    set(gca, 'FontSize', Fsize);
    set(gca,'XtickLabel', [ ] );
    set(gca,'YtickLabel', [ ] );
    if (i == 3 || i == 4)
      Xlab = sprintf('SST: %d K',Sst(i));
      xlabel(Xlab);
    end
    if (i == 1 || i == 3)
      Ylab = sprintf('CCN: %d/cc',Ccn(i));
      ylabel(Ylab);
    end
    % draw the colorbar in the last panel, then output this loop
    % move it to below the panels
    if (i == 4)
      Cbar = colorbar('location', 'SouthOutside');
    end
end

% move the colorbar
% position vector --> [ llx lly width height ]
set(Cbar, 'position', [ 0.1 0.1 0.8 0.05 ]);


saveas(Fig,OutFile);

colormap('default');
close(Fig);


end

