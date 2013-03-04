function [ ] = GenSlopePlots(ConfigFile)
% GenSlopePlots generate slope (POP) plots

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
UndefVal = Config.UndefVal;

Ddir = Config.DiagDir;
Pdir = Config.PlotDir;
% make sure output directory exists
if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

Fsize = 20;

for iplot = 1:length(Config.SlopePlots)
    clear SLOPES;
    clear Y;

    Svar  = Config.SlopePlots(iplot).Svar;
    Cvar  = Config.SlopePlots(iplot).Cvar;
    BLvar = Config.SlopePlots(iplot).BLvar;
    BUvar = Config.SlopePlots(iplot).BUvar;
    Tvar  = Config.SlopePlots(iplot).Tvar;
    Nsamp = Config.SlopePlots(iplot).Nsamp;
    Tcor  = Config.SlopePlots(iplot).Tcor;
    Tmin  = Config.SlopePlots(iplot).Tmin;
    Tmax  = Config.SlopePlots(iplot).Tmax;
    Cmin  = Config.SlopePlots(iplot).Cmin;
    Cmax  = Config.SlopePlots(iplot).Cmax;

    Ptitle = sprintf('%s', Config.SlopePlots(iplot).Title);
    Xlabel = Config.SlopePlots(iplot).Xlabel;
    Ylabel = Config.SlopePlots(iplot).Ylabel;
    OutFile = sprintf('%s/%s', Pdir, Config.SlopePlots(iplot).OutFile);

    fprintf('Generating slope plot:\n');
    fprintf('  Time range for data selection (hr):\n');
    fprintf('    Tmin: %.2f\n', Tmin);
    fprintf('    Tmax: %.2f\n', Tmax);
    fprintf('  Parameters for Fisher z-score screening:\n');
    fprintf('    Number of samples: %d\n', Nsamp);
    fprintf('    Test correlation: %.2f\n', Tcor);
    fprintf('\n');
    
    ips = Config.SlopePlots(iplot).PSnum;
    if (ips == 0)
      fprintf('WARNING: skipping SlopePlot number %d due to no associated PlotSet\n', iplot)
    else
      for icase = 1:Config.PlotSets(ips).Ncases
        Case = Config.PlotSets(ips).Cases(icase).Cname;
        % steal the legend text to generate the Y axis values
        Y(icase) = sscanf(Config.PlotSets(ips).Cases(icase).Legend, '%d');

        % Slopes and correlation coefficients are organized (x,t)
        % If on the first case, also read in lower and upper bin edges, and time
        Hfile = sprintf('%s/%s.h5', Ddir, Case);
        fprintf('  Reading HDF5 file: %s\n', Hfile);
        S = hdf5read(Hfile, Svar);
        C = hdf5read(Hfile, Cvar);
        if (icase == 1)
          % read in bin edge values and time
          BL = hdf5read(Hfile, BLvar);
          BU = hdf5read(Hfile, BUvar);
          T  = hdf5read(Hfile, Tvar)/3600; % hr

          % Create the LWP bin values
          %    Use centers for now
          X = (BL + BU) .* 0.5;

          % find the indices for selecting out a time range
          T1 = find(T >= Tmin, 1, 'first');
          T2 = find(T <= Tmax, 1, 'last');
        end

        % Select out the slope data based on the time range
        S = S(:,T1:T2);
        C = C(:,T1:T2);

        % Select out the slopes that have >= 95% statistical significance
        % Use the fisher z-score to test the null hypothesis that the
        % population correlation is <= Tcor. We want to reject the null
        % hypothesis since we want to keep the slopes from regression fits
        % that indicate that the population correlation > Tcor. For 2-tailed
        % test at 95% significance, need the z-score >= 1.96. Place nans
        % in S where the z-score is < 1.96 so that these slopes do not get
        % figured into the time average.
        S(FisherZscore(Tcor, C, Nsamp) < 1.96) = nan;

        % Save time average of slopes
        % Put the cases in the rows, so that in the contour plot
        % the cases will end up being the y-axis
        SLOPES(icase,:) = nanmean(S,2);
      end
    end


    fprintf('  Writing plot file: %s\n', OutFile);
    Fig = figure;

    contourf(X, Y, SLOPES);
    shading flat;
    set(gca, 'FontSize', Fsize);
    title(Ptitle);
    xlabel(Xlabel);
    ylabel(Ylabel);
    cbar = colorbar;
    set(cbar, 'FontSize', Fsize);
    colormap(redblue);
    caxis([ Cmin Cmax ]);

    saveas(Fig, OutFile);
    close(Fig);
    fprintf('\n');
end
