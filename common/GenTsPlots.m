function [ ] = GenTsPlots(ConfigFile)
% GenTsPlots generate time series plots

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);
    
Pname = Config.Pexp.Ename;
Tstart = Config.Pexp.Tstart;
Tend = Config.Pexp.Tend;

Tdir = Config.TsavgDir;
Pdir = Config.PlotDir;

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
        case 2
            Var = 'horiz_ke';
            Ptitle = sprintf('%s: Total Horizontal Kinetic Energy', Pname);
            Ylim = [ 0 3.5e11 ];
            Ylabel = 'KE (J)';
            OutFile = sprintf('%s/HorizKE.jpg', Pdir);
            LegLoc = 'NorthWest';
        case 3
            Var = 'storm_int';
            Ptitle = sprintf('%s: Storm Intensity Metric', Pname);
            Ylim = [ 0 0.08 ];
            Ylabel = 'Metric';
            OutFile = sprintf('%s/StormInt.jpg', Pdir);
            LegLoc = 'NorthWest';
        otherwise
            fprintf('WARNING: Unrecongnized Ptype');
    end
    
    
    % make sure output directory exists
    if (exist(Pdir, 'dir') ~= 7)
        mkdir(Pdir);
    end
    
    for icase = 1:length(Config.Cases)
        Case = Config.Cases(icase).Cname;

        Hfile = sprintf('%s/%s_%s.h5', Tdir, Var, Case);
        fprintf('Reading HDF5 file: %s\n', Hfile);
        [ Hvar, Rcoords, Zcoords, Tcoords ] = ReadAzavgVar(Hfile, Var);
        TS = squeeze(Hvar);
        
        % smooth with a running mean of length 'Flen'
        [ TsAll(icase,:) ] = SmoothFillTseries(TS, Tlen, Flen);
        
        LegText(icase) = { Config.Cases(icase).Pname };
    end
    
    PlotTseriesSet( Times, TsAll, Ptitle, Ylim, Ylabel, Lcolors, LegText, LegLoc, OutFile );
    fprintf('\n');
end
