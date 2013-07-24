function [] = PlotSampleW250(ConfigFile)

Config = ReadConfig(ConfigFile);

Hdir = 'HDF5';
Pdir = Config.PlotDir;

if (exist(Pdir, 'dir') ~= 7)
  mkdir(Pdir);
end

Tlist = [ 217 289 361 433 ];

% W will be organized as (x,y,t)
Hfile = sprintf('%s/w_z250-z.atex.ccn1600.sst303.gcn10m5.1um-AS-1999-02-10-040000-g1.h5', Hdir);
Hdset = 'w';
fprintf('Reading: %s, Dataset: %s\n', Hfile, Hdset);

W250 = squeeze(hdf5read(Hfile, Hdset));
X = hdf5read(Hfile, 'x_coords');
Y = hdf5read(Hfile, 'y_coords');

Xlabel = 'Longitude';
Ylabel = 'Latitude';
Crange = [ -1.5 1.5 ];

for it = 1:length(Tlist)
  T = Tlist(it);
  Time = ((T-1) * 5) / 60; % T starts at Time = 0, each T is 5 minutes, convert to hours

  Ptitle = sprintf('W: C1600 S303: Z = 250 m, T = %.1f h', Time);
  OutFile = sprintf('%s/Sample_W250_T%d.jpg', Pdir, T);

  Pdata = squeeze(W250(:,:,T))';

  fprintf('Writing: %s\n', OutFile);
  MakeWplotFile(X, Y, Pdata, Ptitle, Xlabel, Ylabel, Crange, OutFile);
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
