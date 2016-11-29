function [ ] = PlotFsFigCfracProfiles()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 13;

  CaseList = {
%    { 'TSD_SAL_DUST'       'SD'   'black' }
%    { 'TSD_NONSAL_DUST'    'NSD'  'red'   }
    { 'TSD_SAL_NODUST'     'SND'  'blue'  }
    { 'TSD_NONSAL_NODUST'  'NSND' 'green' }
    };
  Nc = length(CaseList);

  Zmin =  0; % km
  Zmax = 10; % km

  TstartPsap = 10; % h
  TendPsap   = 30; % h

  % Read in cloud fraction
  for ic = 1:Nc
    Case  = CaseList{ic}{1};
    Label = CaseList{ic}{2};
    Color = CaseList{ic}{3};

    InFname = sprintf('TsAveragedData/hfrac_cloud_%s.h5', Case);
    InSalArFname = sprintf('TsAveragedData/hfrac_sal_ar_cloud_%s.h5', Case);
    InVname  = '/cloud';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    fprintf('Reading: %s (%s)\n', InSalArFname, InVname);

    if (ic == 1)
      X = squeeze(h5read(InFname, '/x_coords')); % km
      Y = squeeze(h5read(InFname, '/y_coords')); % km
      Z = squeeze(h5read(InFname, '/z_coords'))./1000; % km
      T = squeeze(h5read(InFname, '/t_coords'))./3600 - 42; % sim time, h

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % set height indices
      Z1 = find(Z >= Zmin, 1, 'first');
      Z2 = find(Z <= Zmax, 1, 'last');

      % set time indices
      T1 = find(T >= TstartPsap, 1, 'first');
      T2 = find(T <= TendPsap,   1, 'last');

      % set hdf5 hypercube selection
      StartNc = [ 1 1  1  1 ]; % y == 1 is count of cells with cloud
      StartN  = [ 1 2  1  1 ]; % y == 2 is total number of cells
      Count   = [ 1 1 Nz Nt ];

      % create output arrays
      HGHTS = Z(Z1:Z2);
      Nzout = length(HGHTS);

      CFRAC        = zeros([ Nzout Nc ]);
      CFRAC_SAL_AR = zeros([ Nzout Nc ]);
    end

    CloudyCells = squeeze(h5read(InFname, InVname, StartNc, Count));
    TotalCells  = squeeze(h5read(InFname, InVname, StartN,  Count));
    CfracAllT = CloudyCells ./ TotalCells;
    CFRAC(:,ic) = mean(CfracAllT(Z1:Z2,T1:T2),2);

    CloudyCells = squeeze(h5read(InSalArFname, InVname, StartNc, Count));
    TotalCells  = squeeze(h5read(InSalArFname, InVname, StartN,  Count));
    CfracAllT = CloudyCells ./ TotalCells;
    CFRAC_SAL_AR(:,ic)  = mean(CfracAllT(Z1:Z2,T1:T2),2);

    LegText{ic} = Label;
    Colors{ic} = Color;
  end

  fprintf('\n');

  % plot
  Fig = figure;

  Xlim = [ 0 1 ];
  Ylim = [ 0 10 ];

  Paxes = subplot(1,2,1);
  PlotFsFigLine(Paxes, CFRAC, HGHTS, 'a', 'PSAP: G3', 'Cloud Frac.', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'NorthEast', Colors);

  Paxes = subplot(1,2,2);
  PlotFsFigLine(Paxes, CFRAC_SAL_AR, HGHTS, 'b', 'PSAP: SAL\_AR', 'Cloud Frac.', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);


  OutFile = sprintf('%s/FsFigCfracProfiles.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
