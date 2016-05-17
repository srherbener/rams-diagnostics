function [ ] = PlotFsFigCloudRainProfiles()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 13;

  CaseList = {
    { 'TSD_SAL_DUST'       'SD'   'black' }
    { 'TSD_NONSAL_DUST'    'NSD'  'red'   }
    { 'TSD_SAL_NODUST'     'SND'  'blue'  }
    { 'TSD_NONSAL_NODUST'  'NSND' 'green' }
    };
  Nc = length(CaseList);
    

  % Read in dust1, dust2 number concentrations
  for ic = 1:Nc
    Case  = CaseList{ic}{1};
    Label = CaseList{ic}{2};
    Color = CaseList{ic}{3};

    InFname = sprintf('DIAGS/storm_profiles_%s.h5', Case);
    CnumVname  = '/all_rb_s_cloud_num';
    CdiamVname = '/all_rb_s_cloud_diam';
    RmassVname = '/all_rb_s_rain_mass';

    fprintf('Reading: %s (%s)\n', InFname, CnumVname);
    fprintf('Reading: %s (%s)\n', InFname, CdiamVname);
    fprintf('Reading: %s (%s)\n', InFname, RmassVname);
    fprintf('\n');

    if (ic == 1)
      Z = squeeze(h5read(InFname, '/z_coords'))./1000; % km
    end

    CloudNum(:,ic) = squeeze(h5read(InFname, CnumVname));
    CloudDiam(:,ic) = squeeze(h5read(InFname, CdiamVname));
    RainMass(:,ic) = squeeze(h5read(InFname, RmassVname));
    LegText{ic} = Label;
    Colors{ic} = Color;
  end

  % plot
  Fig = figure;

  Ylim = [ 0 5 ];

  Paxes = subplot(1,3,1);
  Xlim = [ 0 80 ];
  PlotFsFigProfile(Paxes, CloudNum, Z, 'a', 'SAL: Rband', 'N_c (cm^-^3)', Xlim, 'Height (km)', Ylim, Fsize, LegText, 'NorthEast', Colors);

  Paxes = subplot(1,3,2);
  Xlim = [ 0 60 ];
  PlotFsFigProfile(Paxes, CloudDiam, Z, 'b', 'SAL: Rband', 'D_c (\mum)', Xlim, 'Height (km)', Ylim, Fsize, LegText, 'none', Colors);

  Paxes = subplot(1,3,3);
  Xlim = [ 0 0.15 ];
  PlotFsFigProfile(Paxes, RainMass, Z, 'c', 'SAL: Rband', 'Rain (g kg^-^1)', Xlim, 'Height (km)', Ylim, Fsize, LegText, 'none', Colors);

  OutFile = sprintf('%s/FsFigCloudRainProfiles.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
