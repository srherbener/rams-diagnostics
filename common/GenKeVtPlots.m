% script to plot KE vs Vt

clear;

[ Config ] = ReadConfig('DiagConfig');

Pname = Config.Pexp.Ename;
Tstart = Config.Pexp.Tstart;
Tend = Config.Pexp.Tend;

Tdir = Config.TsavgDir;
Pdir = Config.PlotDir;
ControlCase = Config.ControlCase;

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

for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;
    Pcase = Config.Cases(icase).Pname;

    KeFile = sprintf('%s/%s_%s.h5', Tdir, KeVar, Case);
    fprintf('Reading HDF5 file: %s\n', KeFile);
    [ KE, Rcoords, Zcoords, Tcoords ] = ReadAzavgVar(KeFile, KeVar);
    [ KeAll(icase,:) ] = SmoothFillTseries(squeeze(KE), Tlen, Flen);
    
    VtFile = sprintf('%s/%s_%s.h5', Tdir, VtVar, Case);
    fprintf('Reading HDF5 file: %s\n', VtFile);
    [ VT, Rcoords, Zcoords, Tcoords ] = ReadAzavgVar(VtFile, VtVar);
    [ VtAll(icase,:) ] = SmoothFillTseries(squeeze(VT), Tlen, Flen);
    
    LegText(icase) = { Pcase };
end

Plot2dSet( VtAll, KeAll, Ptitle, Xlabel, Ylabel, Lcolors, LegText, LegLoc, OutFile );
