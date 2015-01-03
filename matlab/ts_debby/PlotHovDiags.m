function [ ] = PlotHovDiags(ConfigFile)

  [ Config ] = ReadConfig(ConfigFile);

  Pdir = Config.PlotDir;
  Adir = Config.AzavgDir;

  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  VarList = {
    %   FILE        VAR NAME    OUT FILE       Crange     Clevs
    { 'vint_cond' '/vint_cond' 'hov_vint_cond' [ 0  4 ] [ (0:0.2:4) ] }
    { 'pcprate_1' '/pcprate'   'hov_pcprate_1' [ 0 10 ] [ (0:1:10)  ] }
    };


  Nc = length(Config.Cases);
  Nv = length(VarList);

  for icase = 1:Nc
    Case = Config.Cases(icase).Cname;

    fprintf('***********************************************************************\n');
    fprintf('Generating Hovmoller Diagrams\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    for ivar = 1:Nv
      InFprefix  = VarList{ivar}{1};
      InVname    = VarList{ivar}{2};
      OutFprefix = VarList{ivar}{3};
      Crange     = VarList{ivar}{4};
      Clevs      = VarList{ivar}{5};

      InFile = sprintf('%s/%s_%s.h5', Adir, InFprefix, Case);
      fprintf('  Reading: %s (%s)\n', InFile, InVname);

      % Data: (r,t)
      HDATA = squeeze(h5read(InFile, InVname));
      R     = squeeze(h5read(InFile, '/x_coords')) ./ 1000;  % km
      T     = (squeeze(h5read(InFile, '/t_coords')) ./ 3600) - 42.0;  % hr
 
      % Data needs to be transposed for plotting.
      % Replace -999 with nan
      PDATA = HDATA';
      PDATA(PDATA == -999) = nan;

      % Make plot
      OutFile = sprintf('%s/%s_%s.jpg', Pdir, OutFprefix, Case);
      fprintf ('  Writing: %s\n', OutFile);

      Fsize = 22;

      Fig = figure;
      contourf(R, T, PDATA, Clevs, 'LineStyle', 'none');
      set(gca, 'FontSize', Fsize);
      caxis(Crange);
      xlabel('Radius (km)');
      ylabel('Time (h)');
      colorbar;

      saveas(Fig, OutFile);
      close(Fig);

      fprintf('\n');
    end
  end
end
