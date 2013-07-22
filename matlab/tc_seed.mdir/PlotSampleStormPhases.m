function [ ] = PlotSampleStormPhases( ConfigFile )
%PlotSampleStormPhases make plot showing phases of TC
%   This function will read in max Vt time series, plot them out and mark
%   where the rapid intensification (RI) phase and the steady state (SS)
%   phase. The intention of this is to show where these phases are being
%   sampled for diagnostics.

Pset = 'ams2';
Tmin = 24; %hr
Tmax = 90;
Vname = 'max_azwind';
Fprefix = 'max_azwind';
Flen = 5;

Lwidth = 2;
Fsize = 20;
LegLoc = 'SouthEast';
Xlabel = 'Simulation Time (hr)';
Ylabel = 'Vt (m/s)';
Xlim = [ 10 100 ];
Ylim = [ 10 60 ];

% % Rectangular patches for drawing regions where RI and SS exist
% PX = [ 50 90;
%        50 90;
%        70 110;
%        70 110 ];
% PY = [ 25 35;
%        55 55;
%        55 55;
%        25 35 ];
% PC = [ 0 0;
%        0 0;
%        0 0;
%        0 0 ];

% Rectangular patches for drawing regions where RI exist
PX = [ 50;
       50;
       70;
       70 ];
PY = [ 25;
       55;
       55;
       25 ];
PC = [ 0;
       0;
       0;
       0 ];
PG = 0.8; % Gray scale


[ Config ] = ReadConfig(ConfigFile);

Tdir = Config.TsavgDir;
Pdir = Config.PlotDir;

ExpName = Config.ExpName;
Ptitle = sprintf('Maximum Vt');

OutFile = sprintf('%s/ams_StormPhases.jpg', Pdir);

% Make sure PlotDir exists
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

% Find the plot set
ips = find(strcmp(Pset, { Config.PlotSets(:).Name }));

% Read in the Vt and Time data
if (isempty(ips))
    fprintf('ERROR: Cannot find PlotSet "%s" in ConfigFile: %s\n', Pset, ConfigFile);
else
    Ncases = Config.PlotSets(ips).Ncases;
    for icase = 1:Ncases
        Case = Config.PlotSets(ips).Cases(icase).Cname;
        LegText(icase) = { Config.PlotSets(ips).Cases(icase).Legend };
        LineSpecs(icase) = { Config.PlotSets(ips).Cases(icase).Lspec };
        
        Hfile = sprintf('%s/%s_%s.h5', Tdir, Fprefix, Case);
        fprintf('Reading HDF5 file: %s, Dataset: %s\n', Hfile, Vname);
        VT = hdf5read(Hfile, Vname);
        T = hdf5read(Hfile, '/t_coords') / 3600;
        
        % Select data according to Tmin, Tmax
        T1 = find(T >= Tmin, 1, 'first');
        T2 = find(T <= Tmax, 1, 'last');
        
        VT = squeeze(VT(:,:,:,T1:T2));
        T = T(T1:T2);
        
        VTall(icase,:) = SmoothFillTseries(VT, length(VT), Flen);
        if (icase == 1)
            Times = T;
        end
    end
    
    % Construct the plot
    Fig = figure;
    
    % Draw the lines
    plot(Times,VTall(1,:),LineSpecs{1},'LineWidth', Lwidth);
    set(gca, 'FontSize', Fsize);
    set(gca, 'Xlim', Xlim);
    set(gca, 'Ylim', Ylim);
    hold on;
    for iplot = 2:Ncases
        plot(Times,VTall(iplot,:),LineSpecs{iplot}, 'LineWidth', Lwidth);
    end
    
    % plot labels
    legend(LegText, 'Location', LegLoc);
    legend boxoff;
    title(Ptitle);
    xlabel(Xlabel);
    ylabel(Ylabel);
    
    % draw a transparent rectangle for the RI and SS regions
    %R = annotation('rectangle', [.5 0 .5 .5 ] );
    %set(R, 'FaceColor', 'r');
    %set(R, 'FaceAlpha', 0.25);
    %rectangle('Position', [ 50 25 20 30 ], 'LineStyle', '--');
    %rectangle('Position', [ 90 35 20 20 ], 'LineStyle', '--');
    patch(PX, PY, PC, 'FaceColor', [ PG PG PG ], 'LineStyle', '--');
    text(60,25,'RI', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom', 'FontSize', Fsize);
    %text(100,35,'SS', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom', 'FontSize', Fsize);
    
    fprintf('Writing plot file: %s\n', OutFile);
    saveas(Fig, OutFile);
    
    hold off;
    close(Fig);
    
end

end

