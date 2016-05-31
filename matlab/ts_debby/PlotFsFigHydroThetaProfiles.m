function [ ] = PlotFsFigHydroThetaProfiles()

  Pdir = 'Plots';
  if (exist(Pdir, 'dir') ~= 7)
      mkdir(Pdir);
  end
  
  Fsize = 9;

  CaseList = {
    { 'TSD_SAL_DUST'       'SD'   'black' }
    { 'TSD_NONSAL_DUST'    'NSD'  'red'   }
    { 'TSD_SAL_NODUST'     'SND'  'blue'  }
    { 'TSD_NONSAL_NODUST'  'NSND' 'green' }
    };
  Nc = length(CaseList);
    

  % Read in profiles
  for ic = 1:Nc
    Case  = CaseList{ic}{1};
    Label = CaseList{ic}{2};
    Color = CaseList{ic}{3};

    InFname = sprintf('DIAGS/storm_profiles_%s.h5', Case);
    InVname  = '/all_rb_s_cloud_mass';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    CloudMass(:,ic) = squeeze(h5read(InFname, InVname));
    if (ic == 1)
      Z = squeeze(h5read(InFname, '/z_coords'))./1000; % km
    end

    InVname  = '/all_rb_s_rain_mass';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    RainMass(:,ic) = squeeze(h5read(InFname, InVname));

    InVname = '/all_rb_s_pris_mass';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    PrisMass(:,ic) = squeeze(h5read(InFname, InVname));

    InVname = '/all_rb_s_snow_mass';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    SnowMass(:,ic) = squeeze(h5read(InFname, InVname));

    InVname = '/all_rb_s_aggr_mass';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    AggrMass(:,ic) = squeeze(h5read(InFname, InVname));



    InVname = '/all_rb_s_vapor';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    Vapor(:,ic) = squeeze(h5read(InFname, InVname));

    InVname = '/all_rb_s_relhum';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    RelHum(:,ic) = squeeze(h5read(InFname, InVname));

    InVname = '/all_rb_s_ice_sub';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    IceSub(:,ic) = squeeze(h5read(InFname, InVname));

    InVname = '/all_rb_s_liq_evap';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    LiqEvap(:,ic) = squeeze(h5read(InFname, InVname));

    InVname = '/all_rb_s_theta_e';
    fprintf('Reading: %s (%s)\n', InFname, InVname);
    ThetaE(:,ic) = squeeze(h5read(InFname, InVname));

    LegText{ic} = Label;
    Colors{ic} = Color;

    fprintf('\n');
  end

  % plot
  Fig = figure;

  Paxes = subplot(3,2,1);
  Xlim = [ 0 0.3 ];
  Ylim = [ 0 2 ];
  PlotFsFigLine(Paxes, CloudMass, Z, 'a', 'SAP: Rband', 'M_c (g kg^-^1)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'NorthEast', Colors);

  Paxes = subplot(3,2,2);
  Xlim = [ 0 0.1 ];
  Ylim = [ 0 2 ];
  PlotFsFigLine(Paxes, RainMass, Z, 'a', 'SAP: Rband', 'M_r (g kg^-^1)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'NorthEast', Colors);

%  Paxes = subplot(3,2,3);
%  Xlim = [ 0  0.1 ];
%  Ylim = [ 0 15 ];
%  PlotFsFigLine(Paxes, PrisMass, Z, 'b', 'SAP: Rband', 'M_p (g kg^-^1)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);
%
%  Paxes = subplot(3,2,3);
%  Xlim = [ 0  0.2 ];
%  Ylim = [ 0 15 ];
%  PlotFsFigLine(Paxes, AggrMass, Z, 'c', 'SAP: Rband', 'M_a (g kg^-^1)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);
% not much diff in snow mass
%  Xlim = [ 0  0.2 ];
%  Ylim = [ 0 15 ];
%  PlotFsFigLine(Paxes, SnowMass, Z, 'c', 'SAP: Rband', 'M_s (g kg^-^1)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);
%
%  Paxes = subplot(3,2,4);
%  Xlim = [ -1 0 ];
%  Ylim = [ 0 15 ];
%  PlotFsFigLine(Paxes, IceSub, Z, 'd', 'SAP: Rband', 'Ice Sub. (g kg^-^1)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);



  Paxes = subplot(3,2,3);
  Xlim = [ 80 100 ];
  Ylim = [  0   2 ];
  PlotFsFigLine(Paxes, RelHum, Z, 'd', 'SAP: Rband', 'RH (%)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);

  Paxes = subplot(3,2,4);
  Xlim = [ 12 18 ];
  Ylim = [  0  2 ];
  PlotFsFigLine(Paxes, Vapor, Z, 'd', 'SAP: Rband', 'Vapor (g kg^-^1)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);

  Paxes = subplot(3,2,5);
  Xlim = [ -1 0 ];
  Ylim = [ 0 2 ];
  PlotFsFigLine(Paxes, LiqEvap, Z, 'd', 'SAP: Rband', 'Liq. Evap. (g kg^-^1)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);

  Paxes = subplot(3,2,6);
  Xlim = [ 340 355 ];
  Ylim = [ 0 2 ];
  PlotFsFigLine(Paxes, ThetaE, Z, 'd', 'SAP: Rband', '\theta_e (K)', Xlim, 'linear', 1, 'Height (km)', Ylim, 'linear', 1, Fsize, LegText, 'none', Colors);

  OutFile = sprintf('%s/FsFigHydroThetaProfiles.jpg', Pdir);
  fprintf('Writing: %s\n', OutFile);
  saveas(Fig, OutFile);
  close(Fig);
end
