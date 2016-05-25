function [ ] = PlotFsFigDustProfiles()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 13;

  CaseList = {
    { 'TSD_SAL_DUST'       'SD'   'black' }
    { 'TSD_NONSAL_DUST'    'NSD'  'red'   }
%    { 'TSD_SAL_NODUST'     'SND'  'blue'  }
%    { 'TSD_NONSAL_NODUST'  'NSND' 'green' }
    };
  Nc = length(CaseList);
    

  % Read in dust1, dust2 number concentrations
  for ic = 1:Nc
    Case  = CaseList{ic}{1};
    Label = CaseList{ic}{2};
    Color = CaseList{ic}{3};

    InFname = sprintf('DIAGS/storm_profiles_%s.h5', Case);
    D1CoreVname = '/all_core_s_d1_num';
    D2CoreVname = '/all_core_s_d2_num';
    D1RbandVname = '/all_rb_s_d1_num';
    D2RbandVname = '/all_rb_s_d2_num';

    fprintf('Reading: %s (%s)\n', InFname, D1CoreVname);
    fprintf('Reading: %s (%s)\n', InFname, D2CoreVname);
    fprintf('Reading: %s (%s)\n', InFname, D1RbandVname);
    fprintf('Reading: %s (%s)\n', InFname, D2RbandVname);
    fprintf('\n');

    if (ic == 1)
      Z = squeeze(h5read(InFname, '/z_coords'))./1000; % km
    end

    Dust1Core(:,ic) = squeeze(h5read(InFname, D1CoreVname));
    Dust2Core(:,ic) = squeeze(h5read(InFname, D2CoreVname));
    Dust1Rband(:,ic) = squeeze(h5read(InFname, D1RbandVname));
    Dust2Rband(:,ic) = squeeze(h5read(InFname, D2RbandVname));
    LegText{ic} = Label;
    Colors{ic} = Color;
  end

  % plot
  Fig = figure;

  Xlim = [ 0 100 ];
  Ylim = [ 0  15 ];

  Paxes = subplot(2,2,1);
  PlotFsFigLine(Paxes, Dust1Core, Z, 'a', 'SAP: Core', 'Dust1 Num (cm^-^3)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'NorthEast', Colors);

  Paxes = subplot(2,2,2);
  PlotFsFigLine(Paxes, Dust1Rband, Z, 'b', 'SAP: Rband', 'Dust1 Num (cm^-^3)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);

  Paxes = subplot(2,2,3);
  PlotFsFigLine(Paxes, Dust2Core, Z, 'c', 'SAP: Core', 'Dust2 Num (cm^-^3)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);

  Paxes = subplot(2,2,4);
  PlotFsFigLine(Paxes, Dust2Rband, Z, 'd', 'SAP: Rband', 'Dust2 Num (cm^-^3)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);

  OutFile = sprintf('%s/FsFigDustProfiles.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
