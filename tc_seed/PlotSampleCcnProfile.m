function [ ] = PlotSampleCcnProfile(ConfigFile)
% PlotSampleCcnProfile function to plot a sample of the CCN source profile
%

[ Config ] = ReadConfig(ConfigFile);


ControlCase = Config.ControlCase;
AzavgDir = Config.AzavgDir;
PlotDir = Config.PlotDir;

% read in any file in order to get the z coordinate values
Hfile = sprintf('%s/speed_t_%s.h5', AzavgDir, ControlCase);
Z = squeeze(hdf5read(Hfile, '/z_coords'));

% chop off the bottom z value since it is below the surface
Z1 = 2;
Z2 = find(Z <= 10000, 1, 'last');
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
Xlims = [ -10 CCN_VAL*2 ];
Ylims = [ min(Z) max(Z) ];

% Plot
Lwidth = 2;
Fsize = 20;
Pfile = sprintf('%s/SampleCcnProf.jpg', PlotDir);
Ptitle = sprintf('Aerosol Source Profile: %d/cc Example', CCN_VAL);
Xlabel = sprintf('CCN concentration (#/cc)');
Ylabel = sprintf('Height (m)');

Fig = figure;

% data
P = plot(Z,CCN,'LineWidth',Lwidth);
rotate(P, [ 0 0 1 ], 90, [ 0 0 0 ]);
rotate(P, [ 0 1 0 ], 180, [ 0 0 0 ]);
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
