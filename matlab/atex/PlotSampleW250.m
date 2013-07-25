function [] = PlotSampleW250(ConfigFile)

Config = ReadConfig(ConfigFile);

Hdir = 'HDF5';
Pdir = Config.PlotDir;

if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

Tlist = [ 217 289 361 433 ];
Cases = {
    'z.atex.ccn0050.sst298.gcn10m5.1um'
    'z.atex.ccn0400.sst298.gcn10m5.1um'
    'z.atex.ccn1600.sst298.gcn10m5.1um'
    'z.atex.ccn0400.sst293.gcn10m5.1um'
    'z.atex.ccn0400.sst298.gcn10m5.1um'
    'z.atex.ccn0400.sst303.gcn10m5.1um'
    'z.atex.ccn0400.sst298.gcn10m0.3um'
    'z.atex.ccn0400.sst298.gcn10m2.3um'
    'z.atex.ccn0400.sst298.gcn10m4.3um'
    };

CaseNames = {
    'C50 S298'
    'C400 S298'
    'C1600 S298'
    'C400 S293'
    'C400 S298'
    'C400 S303'
    'C400 S298 G0'
    'C400 S298 G2'
    'C400 S298 G4'
    };


Xlabel = 'Longitude';
Ylabel = 'Latitude';
Crange = [ -1.0 1.0 ];

for ic = 1:length(Cases)
    Case = Cases{ic};
    Cname = CaseNames{ic};
    
    % W will be organized as (x,y,t)
    Hfile = sprintf('%s/w_z250-%s-AS-1999-02-10-040000-g1.h5', Hdir, Case);
    Hdset = 'w';
    fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);

    W250 = squeeze(hdf5read(Hfile, Hdset));
    X = hdf5read(Hfile, 'x_coords');
    Y = hdf5read(Hfile, 'y_coords');

    for it = 1:length(Tlist)
      T = Tlist(it);
      Time = ((T-1) * 5) / 60; % T starts at Time = 0, each T is 5 minutes, convert to hours

      Ptitle = sprintf('%s: Z = 250 m, T = %.1f h', Cname, Time);
      OutFile = sprintf('%s/Sample_W250_T%d_%s.jpg', Pdir, T, Case);

      Pdata = squeeze(W250(:,:,T))';

      fprintf('Writing: %s\n', OutFile);
      MakeWplotFile(X, Y, Pdata, Ptitle, Xlabel, Ylabel, Crange, OutFile);
    end
end

end


%%%%%%%%%%%%%%%%%%%%%
function [] = MakeWplotFile(X, Y, Z, Ptitle, Xlabel, Ylabel, Crange, OutFile)

Fig = figure;

PlotW(Fig, X, Y, Z, Ptitle, Xlabel, Ylabel, Crange);

saveas(Fig, OutFile);
close(Fig);

end

%%%%%%%%%%%%%%%%%%%%%%
function [] = PlotW(Fig, X, Y, Z, Ptitle, Xlabel, Ylabel, Crange)

Fsize = 25;

figure(Fig);

contourf(X, Y, Z);
set(gca, 'FontSize', Fsize);
shading flat;
colormap(redblue);
caxis(Crange);
cbar = colorbar;
set(cbar, 'FontSize', Fsize);

title(Ptitle);
xlabel(Xlabel);
ylabel(Ylabel);

end
