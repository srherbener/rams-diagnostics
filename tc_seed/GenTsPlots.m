% script to plot time series of max/min/mean values

clear;

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Cases, Tdirs, Pexp ] = ReadConfig('DiagConfig');
    
Pname = Pexp{1};
Tstart = Pexp{2};
Tend = Pexp{3};

Tdir = 'TsAveragedData';
Pdir = 'plots';

Flen = 5;
Times = (Tstart:Tend);
Tlen = length(Times);

Lcolors = { 'k' 'm' 'b' 'c' 'g' 'y' 'r' };


% Ptype values:
%   1 - max azwind
%   2 - horiz KE
%   3 - storm intensity
for Ptype = 1:3
    switch Ptype
        case 1
            Var = 'max_azwind';
            Ptitle = sprintf('%s: Maximum Azimuthally Averaged Tangential Wind Speed', Pname);
            Ylim = [ 0 80 ];
            Ylabel = 'Wind Speed (m/s)';
            OutFile = sprintf('%s/MaxTanWind.jpg', Pdir);
            LegLoc = 'SouthEast';
            TsType = 'verbatim';
        case 2
            Var = 'horiz_ke';
            Ptitle = sprintf('%s: Total Horizontal Kinetic Energy', Pname);
            Ylim = [ 0 3.5e11 ];
            Ylabel = 'KE (J)';
            OutFile = sprintf('%s/HorizKE.jpg', Pdir);
            LegLoc = 'NorthWest';
            TsType = 'verbatim';
        case 3
            Var = 'storm_int';
            Ptitle = sprintf('%s: Storm Intensity Metric', Pname);
            Ylim = [ 0 0.08 ];
            Ylabel = 'Metric';
            OutFile = sprintf('%s/StormInt.jpg', Pdir);
            LegLoc = 'NorthWest';
            TsType = 'verbatim';
        otherwise
            fprintf('WARNING: Unrecongnized Ptype');
    end
    
    
    % make sure output directory exists
    if (exist(Pdir, 'dir') ~= 7)
        mkdir(Pdir);
    end
    
    for icase = 1:length(Cases)
        Hfile = sprintf('%s/%s_%s.h5', Tdir, Var, Cases{icase});
        fprintf('Reading HDF5 file: %s\n', Hfile);
        [ Hvar, Rcoords, Zcoords, Tcoords ] = ReadAzavgVar(Hfile, Var);
        
        switch TsType
            case 'min'
                TS = squeeze(min(min(Hvar,[],1),[],2))';
            case 'max'
                TS = squeeze(max(max(Hvar,[],1),[],2))';
            case 'mean'
                TS = squeeze(mean(mean(Hvar,1),2))';
            case 'verbatim'
                TS = squeeze(Hvar);
            otherwise
                fprintf('WARNING: Unexpected time series type: %s, skipping plot.\n', TsType);
        end
        
        % smooth with a running mean of length 'Flen'
        [ TsAll(icase,:) ] = SmoothFillTseries(TS, Tlen, Flen);
        
        % Create legend text, the non-control experiment names have the CCN
        % concentration built into their names.
        if (strcmp(Cases{icase},'TCS_CNTL'))
            LegText(icase) = { 'CONTROL' };
        else
            LegText(icase) = { sprintf('CCN: %d/cc', sscanf(Cases{icase},'TCS_GN_C%d')) };
        end
    end
    
    PlotTseriesSet( Times, TsAll, Ptitle, Ylim, Ylabel, Lcolors, LegText, LegLoc, OutFile );
end
