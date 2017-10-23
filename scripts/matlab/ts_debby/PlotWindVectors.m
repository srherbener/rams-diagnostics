function [ ] = PlotWindVectors()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  Cases = {
   'TSD_SAL_DUST'
   'TSD_SAL_NODUST'
   'TSD_NONSAL_DUST'
   'TSD_NONSAL_NODUST'
   };
  Nc = length(Cases);

  CaseNames = {
   'SD'
   'SND'
   'NSD'
   'NSND'
   };

  % Read in sample wind fields from pre sal period (t = 20 hrs) and
  % at z = 3500 m. This is where the enhanced jet appears.
  T1 = 41; % sim time 20 hrs
%  T1 = 61; % sim time 30 hrs
  Z1 = 26; % level nearest 3.5 km (3.494 km)

  LatBounds = [ 7 23 ];
  LonBounds = [ -40 -14 ];

  Inc = 30;
  Fsize = 25;
  Qscale = 1.5;

  for icase = 1:Nc
    Case = Cases{icase};
    CaseName = CaseNames{icase};

    fprintf('*******************************************************\n');
    fprintf('Plotting Wind Vectors for case: %s\n', Case);
    fprintf('\n');

    Ufile = sprintf('HDF5/%s/HDF5/u-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    Uvname = '/u';
    Vfile = sprintf('HDF5/%s/HDF5/v-%s-AS-2006-08-20-120000-g3.h5', Case, Case);
    Vvname = '/v';

    fprintf('  Reading: %s (%s)\n', Ufile, Uvname);
    fprintf('  Reading: %s (%s)\n', Vfile, Vvname);
    fprintf('\n');

    % Read in vars
    U_DS = ncgeodataset(Ufile);
    V_DS = ncgeodataset(Vfile);

    U_VAR = U_DS.geovariable(Uvname);
    V_VAR = V_DS.geovariable(Vvname);

    X_VAR = U_DS.geovariable('/x_coords');
    Y_VAR = U_DS.geovariable('/y_coords');
    Z_VAR = U_DS.geovariable('/z_coords');
    T_VAR = U_DS.geovariable('/t_coords');

    X = X_VAR.data(:);
    Y = Y_VAR.data(:);
    U = squeeze(U_VAR.data(T1,Z1,:,:));
    V = squeeze(V_VAR.data(T1,Z1,:,:));

    Nx = length(X);
    Ny = length(Y);

    % U and V are (y,x) now which is what is needed for the vector plot.
    % U and V are too dense, so take every nth point.
    X = X(1:Inc:Nx);
    Y = Y(1:Inc:Ny);
    U = U(1:Inc:Ny,1:Inc:Nx);
    V = V(1:Inc:Ny,1:Inc:Nx);

    % create a mesh for m_quiver
    [ MX MY ] = meshgrid(X, Y);

    % Make a vector plot using quiver command
    OutFile = sprintf('%s/FigWindVectors_%s.jpg', Pdir, Case);
    fprintf('  Writing: %s\n', OutFile);
    fprintf('\n');

    Fig = figure;
    set(gca, 'FontSize', Fsize);

    hold on;
    m_proj('miller', 'lat', LatBounds, 'long', LonBounds);
    m_coast('color', 'k'); % k --> black
    m_grid('linestyle','none','box','fancy','tickdir','out');
    m_quiver(MX, MY, U, V, Qscale);
    hold off;

    saveas(Fig,OutFile);
    close(Fig);
  end
end
