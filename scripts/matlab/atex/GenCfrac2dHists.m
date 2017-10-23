function [ ] = GenCfrac2dHists(ConfigFile)
% GenCfrac2dHists generate time series of 2D histograms involving cloud fraction

% Read the config file to get the structure of how the data is laid out in
% the file system.
[ Config ] = ReadConfig(ConfigFile);

Tdir = Config.TsavgDir;
Ddir = Config.DiagDir;

if(exist(Ddir,'dir') ~= 7)
  mkdir(Ddir);
end

CfracFprefix = 'hda_cloud_frac';
CfracVar = 'hda_cloud_mask';
OutFprefix = 'hist2d_cfrac';

VarFprefixes = {
  'hist_cdepth'
  'hist_ctoptempc'
  };

InVars = {
  'hist_cloud_depth'
  'hist_ctop_tempc'
  };

OutVars = {
  'cdepth'
  'ctoptempc'
  };

CfracBins = 0:0.05:1.0;

% Selection for final histograms
Xselect = {
  '0:4000:1'
  '-5:20:1'
  };

CfracSelect = '0:1:1';

Tmin = 12;
Tmax = 36;

% Find and replace underscores in Ptitle, Ylabel with blank spaces
for icase = 1:length(Config.Cases)
    Case = Config.Cases(icase).Cname;

    for ivar = 1: length(InVars)
      InVar = InVars{ivar};
      Fprefix = VarFprefixes{ivar};
      OutVar = OutVars{ivar};
      Temp = sscanf(Xselect{ivar}, '%f:%f:%f', 3);
      Xmin   = Temp(1);
      Xmax   = Temp(2);
      Xgroup = Temp(3);
      Temp = sscanf(CfracSelect, '%f:%f:%f', 3);
      Ymin   = Temp(1);
      Ymax   = Temp(2);
      Ygroup = Temp(3);

      CfracFile = sprintf('%s/%s_%s.h5', Tdir, CfracFprefix, Case);
      InFile = sprintf('%s/%s_%s.h5', Tdir, Fprefix, Case);
      OutFile = sprintf('%s/%s_%s_%s.h5', Ddir, OutFprefix, OutVar, Case);
  
      fprintf('***************************************************************\n');
      fprintf('Generating 2D histogram:\n');
      fprintf('  Case: %s\n', Case);
      fprintf('  Input cloud fraction file: %s\n', CfracFile);
      fprintf('    Var name: %s\n', CfracVar);
      fprintf('  Input var file: %s\n', InFile);
      fprintf('    Var name: %s\n', InVar);
      fprintf('  Data selection:\n');
      fprintf('    Xrange, Xgroup: [ %.2f %.2f ], %d\n', Xmin, Xmax, Xgroup);
      fprintf('    Yrange, Ygroup: [ %.2f %.2f ], %d\n', Ymin, Ymax, Ygroup);
      fprintf('    Trange: [ %.2f %.2f ]\n', Tmin, Tmax);
      fprintf('\n');
  
      % Read in the the cloud fraction and the other variable data.
      %   CFRAC will be (Nt)
      %   VAR will be (Nb,Nt)
      CFRAC  = squeeze(hdf5read(CfracFile, CfracVar));
      VAR    = squeeze(hdf5read(InFile, InVar));

      % Grab coordidnates from the input variable, except for Y which will be set to the cloud
      % fraction bins.
      X = hdf5read(InFile, 'x_coords');
      Y = CfracBins;
      T = hdf5read(InFile, 't_coords')/3600; % hr

      Nx = length(X);
      Ny = length(Y);
      Nt = length(T);

      % Create the 2d histogram. For each time step, bin the cloud fraction, then add the
      % counts in from the input variable in the appropriate y location. GenCountBins wants
      % the input to be (x,y,z,t) so add in a dummy z dimention.
      C = zeros(Nx, Ny, 1, Nt);
      for i = 1:length(T)
        ybin = find(CFRAC(i) >= Y, 1, 'last');
        C(:,ybin,1,i) = C(:,ybin,1,i) + VAR(:,i);
      end

      % Do the selection and bin combination
      T1 = find(T >= Tmin, 1, 'first');
      T2 = find(T <= Tmax, 1, 'last');

      [ COUNTS, XL, XU, YL, YU ] = GenCountBins(C(:,:,:,T1:T2), X, Y, Xmin, Xmax, Xgroup, Ymin, Ymax, Ygroup);

      % If the time range was more than a single point, then sum up counts across the time
      % dimension. Time will be the last dimension.
      if ((T2-T1) > 0)
        COUNTS = sum(COUNTS, ndims(COUNTS));
      end

      % output
      fprintf('Writing: %s\n', OutFile);
      hdf5write(OutFile, '/COUNTS', COUNTS);
      hdf5write(OutFile, 'XL', XL, 'WriteMode', 'append');
      hdf5write(OutFile, 'XU', XU, 'WriteMode', 'append');
      hdf5write(OutFile, 'YL', YL, 'WriteMode', 'append');
      hdf5write(OutFile, 'YU', YU, 'WriteMode', 'append');
      fprintf('\n');
  end

end
