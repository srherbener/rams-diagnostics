% script to plot time series of max/min/mean values

clear;

% set this to select different plots
%   1 - max azwind
%   2 - horiz KE
%   3 - storm intensity
Ptype = 1;

Tdir = 'TsAveragedData';
Pdir = 'plots';

Flen = 5;

Times = (24:144);
Tlen = length(Times);

Lcolors = { 'k' 'm' 'b' 'c' 'g' 'y' 'r' };


switch Ptype
    case 1
        Var = 'max_azwind';
        Ptitle = 'Maximum Azimuthally Averaged Tangential Wind Speed';
        Ylim = [ 0 80 ];
        Ylabel = 'Wind Speed (m/s)';
        OutFile = sprintf('%s/MaxTanWind.jpg', Pdir);
        LegLoc = 'SouthEast';
        TsType = 'verbatim';
    case 2
        Var = 'horiz_ke';
        Ptitle = 'Total Horizontal Kinetic Energy';
        Ylim = [ 0 3.5e11 ];
        Ylabel = 'KE (J)';
        OutFile = sprintf('%s/HorizKE.jpg', Pdir);
        LegLoc = 'NorthWest';
        TsType = 'verbatim';
    case 3
        Var = 'storm_int';
        Ptitle = 'Storm Intensity Metric';
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

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Exps, Tdirs ] = ReadConfig('Config');

for iexp = 1:length(Exps)
    Hfile = sprintf('%s/%s_%s.h5', Tdir, Var, char(Exps(iexp)));
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
    [ TsAll(iexp,:) ] = SmoothFillTseries(TS, Tlen, Flen);
    
    % Create legend text, the non-control experiment names have the CCN
    % concentration built into their names.
    if (strcmp(char(Exps(iexp)),'TCS_CNTL'))
        LegText(iexp) = { 'CONTROL' };
    else
        LegText(iexp) = { sprintf('CCN: %d/cc', sscanf(char(Exps(iexp)),'TCS_GN_C%d')) };
    end
end

PlotTseriesSet( Times, TsAll, Ptitle, Ylim, Ylabel, Lcolors, LegText, LegLoc, OutFile );
