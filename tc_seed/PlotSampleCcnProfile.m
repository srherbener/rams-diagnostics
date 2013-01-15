function [ ] = PlotSampleCcnProfile(ConfigFile)
% PlotSampleCcnProfile function to plot a sample of the CCN source profile
%

[ Config ] = ReadConfig(ConfigFile);


ControlCase = Config.ControlCase;
AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;

% read in any file in order to get the z coordinate values
Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, ControlCase);
Z = squeeze(hdf5read(Hfile, '/z_coords')) / 1000; % km

% chop off the bottom z value since it is below the surface
Z1 = 2;
Z2 = find(Z <= 10, 1, 'last');
Z = Z(Z1:Z2);
Nz = length(Z);

% in RAMSIN, CCN_H is 27 and CCN_TAPER is 31 but we've chopped
% off the first Z value which means we need to decrement CCN_H and
% CCN_TAPER by one.
CCN_VAL = 500;
CCN_H = 26;
CCN_TAPER = 30;

CCN = zeros(1,Nz);
for i = 1:Nz
    if (i <= CCN_H)
        CCN(i) = CCN_VAL;
    else
        if (i <= CCN_TAPER)
            CCN(i) = ((CCN_TAPER - i) / 4) * CCN_VAL;
        else
            CCN(i) = 0;
        end
    end
end

% pick x range so that CCN_VAL is in the middle
Xlims = [ -50 CCN_VAL*2 ];
Ylims = [ min(Z) max(Z) ];

% Plot
Lwidth = 2;
Fsize = 20;
Pfile = sprintf('%s/SampleCcnProf.fig', PlotDir);
Ptitle = sprintf('Aerosol Source Profile: %d/cc Example', CCN_VAL);
Xlabel = sprintf('Aerosol concentration (#/cc)');
Ylabel = sprintf('Height (km)');

Fig = figure;

% data
plot(CCN,Z,'LineWidth',Lwidth, 'Color', 'k');
xlim(Xlims);
ylim(Ylims);

% axes
set (gca, 'FontSize', Fsize);
title(Ptitle);
xlabel(Xlabel);
ylabel(Ylabel);

saveas(Fig, Pfile);
close(Fig);
end
