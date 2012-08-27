% script to plot KE vs Vt

clear;

[ Cases, Tdirs, Pexp ] = ReadConfig('DiagConfig');

Pname = Pexp{1};
Tstart = Pexp{2};
Tend = Pexp{3};

Tdir = 'TsAveragedData';
Pdir = 'plots';

KeVar = 'horiz_ke';
VtVar = 'max_azwind';

Ptitle = sprintf('%: Kinetic Energy vs Maximum Vt', Pname);
Xlabel = 'Wind Speed (m/s)';
Ylabel = 'KE (J)';
OutFile = sprintf('%s/KeVt.jpg', Pdir);
LegLoc = 'NorthWest';

Flen = 5;

Times = (Tstart:Tend);
Tlen = length(Times);

Lcolors = { 'k' 'm' 'b' 'c' 'g' 'y' 'r' };


% make sure output directory exists
if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
end

for icase = 1:length(Cases)
    KeFile = sprintf('%s/%s_%s.h5', Tdir, KeVar, Cases{icase});
    fprintf('Reading HDF5 file: %s\n', KeFile);
    [ KE, Rcoords, Zcoords, Tcoords ] = ReadAzavgVar(KeFile, KeVar);
    [ KeAll(icase,:) ] = SmoothFillTseries(squeeze(KE), Tlen, Flen);
    
    VtFile = sprintf('%s/%s_%s.h5', Tdir, VtVar, Cases{icase});
    fprintf('Reading HDF5 file: %s\n', VtFile);
    [ VT, Rcoords, Zcoords, Tcoords ] = ReadAzavgVar(VtFile, VtVar);
    [ VtAll(icase,:) ] = SmoothFillTseries(squeeze(VT), Tlen, Flen);
    
    % Create legend text, the non-control case names have the CCN
    % concentration built into their names.
    if (strcmp(Cases{icase},'TCS_CNTL'))
        LegText(icase) = { 'CONTROL' };
    else
        LegText(icase) = { sprintf('CCN: %d/cc', sscanf(Cases{icase},'TCS_GN_C%d')) };
    end
end

Plot2dSet( VtAll, KeAll, Ptitle, Xlabel, Ylabel, Lcolors, LegText, LegLoc, OutFile );
