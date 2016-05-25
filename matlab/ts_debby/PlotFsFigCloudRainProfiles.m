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
    InVname  = '/all_rb_s_cloud_num';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    CloudNum(:,ic) = squeeze(h5read(InFname, InVname));
    if (ic == 1)
      Z = squeeze(h5read(InFname, '/z_coords'))./1000; % km
    end

    InVname = '/all_rb_s_cloud_diam';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    CloudDiam(:,ic) = squeeze(h5read(InFname, InVname));

    InFname = sprintf('DIAGS/hist_meas_az_pcprate_%s.h5', Case);
    InVname = '/hist_all_rb_s_pcprate';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    PrecipRate(:,ic) = squeeze(h5read(InFname, InVname));
    if (ic == 1)
      BINS = squeeze(h5read(InFname, '/y_coords'));
    end

    LegText{ic} = Label;
    Colors{ic} = Color;

    fprintf('\n');
  end

  % Normalize the count to the maximum count in the precip rate
  PrecipRate = PrecipRate ./ max(PrecipRate(:));

  % plot
  Fig = figure;

  Paxes = subplot(2,2,1);
  Xlim = [ 0 30 ];
  Ylim = [ 0 5 ];
  PlotFsFigLine(Paxes, CloudNum, Z, 'a', 'SAP: Rband', 'N_c (cm^-^3)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'NorthEast', Colors);

  Paxes = subplot(2,2,2);
  Xlim = [ 30 50 ];
  PlotFsFigLine(Paxes, CloudDiam, Z, 'b', 'SAP: Rband', 'D_c (\mum)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 0, Fsize, LegText, 'none', Colors);

  Paxes = subplot(2,2,[3 4]);
  Xlim = [ 1e-4 1e2 ];
  Ylim = [ 1e-4   1 ];
  PlotFsFigLine(Paxes, BINS, PrecipRate, 'c', 'SAP: Rband', 'Precip. Rate (mm h^-^1)', Xlim, 'log', 1, 'Norm. Count', Ylim, 'log', 1, Fsize, LegText, 'none', Colors);

  OutFile = sprintf('%s/FsFigCloudRainProfiles.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
