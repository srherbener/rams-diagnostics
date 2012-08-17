% script to plot time series of max/min/mean values

clear;

% time series calculation specs
Var = 'speed_t';
Adir = 'AzAveragedData';
TsType = 'max';

% plotting specs
Times = (24:96);
Ptitle = 'Maximum Azimuthally Averaged Tangential Wind Speed';
Ylim = [ 0 80 ];
Ylabel = 'Wind Speed (m/s)';
Lcolors = { 'k' 'm' 'b' 'c' 'g' 'y' 'r' };
%LegText = { 'CCN: 50/cc' 'CCN: 100/cc' 'CCN: 150/cc' 'CCN: 200/cc' 'CCN: 500/cc' 'CCN: 1000/cc' 'CCN: 2000/cc' };
LegLoc = 'SouthEast';
OutFile = 'PLOTS/MaxTanWind.jpg';


% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Exps, Tdirs ] = ReadConfig('Config');


for iexp = 1:length(Exps)
    Hfile = sprintf('%s/%s_%s.h5', Adir, Var, char(Exps(iexp)));
    fprintf('Reading HDF5 file: %s\n', Hfile);
    [ Hvar, Rcoords, Zcoords, Tcoords ] = ReadAzavgVar(Hfile, Var);
    
    switch TsType
        case 'min'
            TS = squeeze(min(min(Hvar,[],1),[],2))';
        case 'max'
            TS = squeeze(max(max(Hvar,[],1),[],2))';
        case 'mean'
            TS = squeeze(mean(mean(Hvar,1),2))';
        otherwise
            fprintf('WARNING: Unexpected time series type: %s, skipping plot.\n', TsType);
    end
    
    %if (iexp == 1)
    %    TS_CNTL = TS;
    %else
    %  TsAll(iexp-1,:) = TS - TS_CNTL;
    %end
    TsAll(iexp,:) = TS;
    
    % Create legend text, the non-control experiment names have the CCN
    % concentration built into their names.
    if (strcmp(char(Exps(iexp)),'TCS_CNTL'))
        LegText(iexp) = { 'CONTROL' };
    else
        LegText(iexp) = { sprintf('CCN: %d/cc', sscanf(char(Exps(iexp)),'TCS_GN_C%d')) };
    end
end

% Smooth the data: use a running mean of length 5
Flen = 5;
for iexp = 1:size(TsAll,1)
    TsSmooth(iexp,:) = filtfilt(ones(1,Flen)/Flen,1,double(TsAll(iexp,:)));
end


PlotTseriesSet2( Times, TsSmooth, Ptitle, Ylim, Ylabel, Lcolors, LegText, LegLoc, OutFile );
