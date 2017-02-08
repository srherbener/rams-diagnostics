function [ ] = PlotFsFigCloudRainProfiles()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 9;

  CaseList = {
    { 'TSD_SAL_DUST'       'SD'   'black'      }
    { 'TSD_NONSAL_DUST'    'NSD'  'magenta'    }
    { 'TSD_SAL_NODUST'     'SND'  'cyan'       }
    { 'TSD_NONSAL_NODUST'  'NSND' 'sandybrown' }
    };
  Nc = length(CaseList);
    

  % Read in dust1, dust2 number concentrations
  for ic = 1:Nc
    Case  = CaseList{ic}{1};
    Label = CaseList{ic}{2};
    Color = CaseList{ic}{3};

    InFname = sprintf('DIAGS/storm_profiles_%s.h5', Case);
    InVname  = '/all_rb_s_cloud_num';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    CloudNum(:,ic) = squeeze(h5read(InFname, InVname));
    if (ic == 1)
      Z = squeeze(h5read(InFname, '/z_coords'))./1000; % km
    end

    InVname = '/all_rb_s_cloud_diam';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    CloudDiam(:,ic) = squeeze(h5read(InFname, InVname));

    InVname = '/all_rb_s_cloud_evap';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    CloudEvap(:,ic) = squeeze(h5read(InFname, InVname));

    InVname  = '/all_rb_s_rain_num';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    RainNum(:,ic) = squeeze(h5read(InFname, InVname));

    InVname = '/all_rb_s_rain_diam';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    RainDiam(:,ic) = squeeze(h5read(InFname, InVname));

    InVname = '/all_rb_s_rain_evap';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    RainEvap(:,ic) = squeeze(h5read(InFname, InVname));

    LegText{ic} = Label;
    Colors{ic} = Color;

    fprintf('\n');
  end

  % plot
  Fig = figure;

  Paxes = subplot(3,2,1);
  Xlim = [ 0 30 ];
  Ylim = [ 0 5 ];
  PlotFsFigLine(Paxes, CloudNum, Z, 'a', 'SAP: Rband', 'N_c (cm^-^3)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'NorthEast', Colors);

  Paxes = subplot(3,2,3);
  Xlim = [ 30 50 ];
  PlotFsFigLine(Paxes, CloudDiam, Z, 'c', '', 'D_c (\mum)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);

  Paxes = subplot(3,2,5);
  Xlim = [ -2.5 0 ];
  PlotFsFigLine(Paxes, CloudEvap, Z, 'e', '', 'Cloud Evap. (g kg^-^1)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);


  Paxes = subplot(3,2,2);
  Xlim = [ 0 0.05 ];
  PlotFsFigLine(Paxes, RainNum, Z, 'b', '', 'N_r (cm^-^3)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);

  Paxes = subplot(3,2,4);
  Xlim = [ 0.1 0.3 ];
  PlotFsFigLine(Paxes, RainDiam, Z, 'd', '', 'D_r (mm)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);

  Paxes = subplot(3,2,6);
  Xlim = [ -1 0 ];
  PlotFsFigLine(Paxes, RainEvap, Z, 'f', '', 'Rain Evap. (g kg^-^1)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);


  OutFile = sprintf('%s/FsFigCloudRainProfiles.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
