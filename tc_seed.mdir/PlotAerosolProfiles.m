function [ ] = PlotAerosolProfiles(ConfigFile)
% PlotCcnProfiles function to plot CCN profiles used the aerosol source
%

[ Config ] = ReadConfig(ConfigFile);


ControlCase = Config.ControlCase;
AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;

DoingPresentation = false; % make figure for paper
%DoingPresentation = true; % make figure for presentation

% read in any file in order to get the z coordinate values
Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, ControlCase);
Z = squeeze(hdf5read(Hfile, '/z_coords')) / 1000; % km

% replace the bottom z value with zero since it is below the surface
Z(1) = 0;
Z1 = 1;
Z2 = find(Z <= 10, 1, 'last');
Z = Z(Z1:Z2);
Nz = length(Z);

% in RAMSIN, CCN_H is 27 and CCN_TAPER is 31
CCN_VAL = [ 100 500 1000 2000 ];
CCN_H = 27;
CCN_TAPER = 31;
Np = length(CCN_VAL);

% the profile actually starts with Z(2), but this is 25m which on the
% scale of 10km vertical is essentially zero, so just
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
if (DoingPresentation)
  Lwidth = 2;
  Fsize = 25;
  Pfile = sprintf('%s/TCWS0513_AerosolProfiles.jpg', PlotDir);
else
  Lwidth = 3;
  Fsize = 45;
  Pfile = sprintf('%s/AerosolProfiles.jpg', PlotDir);
end

Xlabel = sprintf('Concentration (#/cc)');
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
set(gca, 'LineWidth', 2);
set(gca, 'TickLength', [ 0.015 0.015 ]);

if (DoingPresentation)
  title('Aerosol Source Profiles');
else
  % The title is in a box that adjusts to the amount of characters in
  % the title. Ie, it doesn't do any good to do Left/Center/Right
  % alignment. But, the entire box can be moved to the left side of the
  % plot.
  T = title('(b)');
  set(T, 'Units', 'Normalized');
  set(T, 'HorizontalAlignment', 'Left');
  Tpos = get(T, 'Position');
  Tpos(1) = 0; % line up with left edge of plot area
  set(T, 'Position', Tpos);
end

xlabel(Xlabel);
ylabel(Ylabel);
legend(Ltext, 'Location', 'NorthEast', 'FontSize', 25);
legend boxoff;

if (~DoingPresentation)
% Fix up the positioning
  Ppos = get(gca, 'Position'); % position of plot area
  Ppos(1) = Ppos(1) * 1.00;
  Ppos(2) = Ppos(2) * 0.90;
  Ppos(3) = Ppos(3) * 0.90;
  Ppos(4) = Ppos(4) * 0.90;
  set(gca, 'Position', Ppos);
end

saveas(Fig, Pfile);
close(Fig);
end
