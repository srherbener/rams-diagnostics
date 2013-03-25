function [ ] = PlotCcnProfiles(ConfigFile)
% PlotCcnProfiles function to plot CCN profiles used the aerosol source
%

[ Config ] = ReadConfig(ConfigFile);


ControlCase = Config.ControlCase;
AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;

% read in any file in order to get the z coordinate values
Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, ControlCase);
Z = squeeze(hdf5read(Hfile, '/z_coords')) / 1000; % km

% replace the bottom z value with zero since it is below the surface
Z(1) = 0;
Z1 = 1;
Z2 = find(Z <= 8, 1, 'last');
Z = Z(Z1:Z2);
Nz = length(Z);

% in RAMSIN, CCN_H is 27 and CCN_TAPER is 31
CCN_VAL = [ 100 500 1000 2000 ];
CCN_H = 27;
CCN_TAPER = 31;
Np = length(CCN_VAL);

% the profile actually starts with Z(2), but this is 25m which on the
% scale of 8km vertical is essentially zero, so just
% let the profile go to the bottom
CCN = zeros(Np,Nz);
for i = 1:Nz
    if (i <= CCN_H)
        CCN(:,i) = CCN_VAL;
    else
        if (i <= CCN_TAPER)
            CCN(:,i) = ((CCN_TAPER - i) ./ 4) .* CCN_VAL;
        else
            CCN(:,i) = 0;
        end
    end
end

% pick x range so that CCN_VAL is in the middle
Xlims = [ -50 2500 ];
Ylims = [ min(Z) max(Z) ];

Lstyles = { '-' '--' ':' '-.' };
Colors = { 'k' 'k' 'k' 'k' };
Ltext = { 'AS0100' 'AS0500' 'AS1000' 'AS2000' };

% Plot
Lwidth = 2;
Fsize = 25;
Pfile = sprintf('%s/AerosolProfiles.jpg', PlotDir);
Ptitle = sprintf('Aerosol Source Profile: %d/cc Example', CCN_VAL);
Xlabel = sprintf('Aerosol concentration (#/cc)');
Ylabel = sprintf('Height (km)');

Fig = figure;

% data
for i = 1:Np
  plot(CCN(i,:), Z, 'LineWidth', Lwidth, 'Color', Colors{i}, 'LineStyle', Lstyles{i});
  if (i == 1)
    hold on;
  end
end
xlim(Xlims);
ylim(Ylims);

% axes
set (gca, 'FontSize', Fsize);
%title(Ptitle);
xlabel(Xlabel);
ylabel(Ylabel);
legend(Ltext, 'Location', 'NorthEast');
legend boxoff;


saveas(Fig, Pfile);
close(Fig);
end
