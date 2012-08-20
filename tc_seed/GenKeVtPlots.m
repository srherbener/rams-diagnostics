% script to plot KE vs Vt

clear;


Tdir = 'TsAveragedData';
Pdir = 'plots';

KeVar = 'horiz_ke';
VtVar = 'max_azwind';

Ptitle = 'Kinetic Energy vs Maximum Vt';
Xlabel = 'Wind Speed (m/s)';
Ylabel = 'KE (J)';
OutFile = sprintf('%s/KeVt.jpg', Pdir);
LegLoc = 'NorthWest';

Flen = 5;

Times = (24:144);
Tlen = length(Times);

Lcolors = { 'k' 'm' 'b' 'c' 'g' 'y' 'r' };


% make sure output directory exists
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Exps, Tdirs ] = ReadConfig('Config');

for iexp = 1:length(Exps)
    KeFile = sprintf('%s/%s_%s.h5', Tdir, KeVar, char(Exps(iexp)));
    fprintf('Reading HDF5 file: %s\n', KeFile);
    [ KE, Rcoords, Zcoords, Tcoords ] = ReadAzavgVar(KeFile, KeVar);
    [ KeAll(iexp,:) ] = SmoothFillTseries(squeeze(KE), Tlen, Flen);
    
    VtFile = sprintf('%s/%s_%s.h5', Tdir, VtVar, char(Exps(iexp)));
    fprintf('Reading HDF5 file: %s\n', VtFile);
    [ VT, Rcoords, Zcoords, Tcoords ] = ReadAzavgVar(VtFile, VtVar);
    [ VtAll(iexp,:) ] = SmoothFillTseries(squeeze(VT), Tlen, Flen);
    
    % Create legend text, the non-control experiment names have the CCN
    % concentration built into their names.
    if (strcmp(char(Exps(iexp)),'TCS_CNTL'))
        LegText(iexp) = { 'CONTROL' };
    else
        LegText(iexp) = { sprintf('CCN: %d/cc', sscanf(char(Exps(iexp)),'TCS_GN_C%d')) };
    end
end

Plot2dSet( VtAll, KeAll, Ptitle, Xlabel, Ylabel, Lcolors, LegText, LegLoc, OutFile );
