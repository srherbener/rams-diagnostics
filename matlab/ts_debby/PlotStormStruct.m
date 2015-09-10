function [ ] = PlotStormStruct()

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
   'SAL\_DUST'
   'SAL\_NODUST'
   'NONSAL\_DUST'
   'NONSAL\_NODUST'
   };

  % Two time periods:
  %    Pre SAL: T = 10 to 30 h (after RI, before encounter SAL)
  %        SAL: T = 40 to 60 h (during SAL)

  PreSalTstart = 10;
  PreSalTend   = 30;
  SalTstart    = 40;
  SalTend      = 60;

  % Trim down spatial range to focus on just the storm
  Zbot = 0;  % km
  Ztop = 10; % km

  Rinside  =   0; % km
  Routside = 250; % km

  Fsize = 22;

  VtClim = [ 0 15 ];
  UpClevs = [ 0:0.02:2 ];
  DnClevs = [ -2:0.02:0 ];

  % Want azimuthally averaged tangential speed and vertical velocity for each panel
  VtFileTemplate  = 'DIAGS/hist_meas_speed_<CASE>.h5';
  VtVname = '/avg_speed_t_wm';

  UpFileTemplate  = 'DIAGS/hist_meas_w_<CASE>.h5';
  UpVname = '/avg_updraft_wm';

  DnFileTemplate  = 'DIAGS/hist_meas_w_<CASE>.h5';
  DnVname = '/avg_dndraft_wm';

  PreSalOfileTemplate = 'Plots/PreSalStormStruct_<CASE>.jpg';
  SalOfileTemplate = 'Plots/SalStormStruct_<CASE>.jpg';

  fprintf('Generating storm structure plots\n');
  fprintf('\n');
  for ic = 1:Nc
    Case = Cases{ic};
    CaseName = CaseNames{ic};

    fprintf('**********************************************************\n');
    fprintf('Case: %s\n', Case);
    fprintf('\n');

    VtFile  = regexprep(VtFileTemplate, '<CASE>', Case);
    UpFile  = regexprep(UpFileTemplate, '<CASE>', Case);
    DnFile  = regexprep(DnFileTemplate, '<CASE>', Case);

    PreSalOutFile = regexprep(PreSalOfileTemplate, '<CASE>', Case);
    SalOutFile    = regexprep(SalOfileTemplate, '<CASE>', Case);

    fprintf('  Reading: %s (%s)\n', VtFile, VtVname);
    fprintf('  Reading: %s (%s)\n', UpFile, UpVname);
    fprintf('\n');

    % Read in vars. Vars are organized as (r,z,t) where
    % r is the radius, z is the height, and t is time.
    VT = squeeze(h5read(VtFile, VtVname));
    R  = squeeze(h5read(VtFile, '/x_coords')) ./ 1000; % radius in km
    Z  = squeeze(h5read(VtFile, '/z_coords')) ./ 1000; % height in km
    T  = (squeeze(h5read(VtFile, '/t_coords')) ./ 3600) - 42; % time in hours, 0h -> sim time 42h
    UP = squeeze(h5read(UpFile, UpVname));
    DN = squeeze(h5read(DnFile, DnVname));

    % Find indicies corresponding to selection
    R1 = find(R >= Rinside,  1, 'first');
    R2 = find(R <= Routside, 1, 'last');

    Z1 = find(Z >= Zbot, 1, 'first');
    Z2 = find(Z <= Ztop, 1, 'last');

    PT1 = find(T >= PreSalTstart, 1, 'first');
    PT2 = find(T <= PreSalTend,   1, 'last');
    
    ST1 = find(T >= SalTstart, 1, 'first');
    ST2 = find(T <= SalTend,   1, 'last');

    % vars will be (r,z) after mean is performed
    P_VT = squeeze(mean(VT(R1:R2,Z1:Z2,PT1:PT2),3));
    P_UP = squeeze(mean(UP(R1:R2,Z1:Z2,PT1:PT2),3));
    P_DN = squeeze(mean(DN(R1:R2,Z1:Z2,PT1:PT2),3));

    S_VT = squeeze(mean(VT(R1:R2,Z1:Z2,ST1:ST2),3));
    S_UP = squeeze(mean(UP(R1:R2,Z1:Z2,ST1:ST2),3));
    S_DN = squeeze(mean(DN(R1:R2,Z1:Z2,ST1:ST2),3));

    % adjust the coordinates to match what was selected out of the vars
    R = R(R1:R2);
    Z = Z(Z1:Z2);


    % plot the pre SAL data
    fprintf('  Writing: %s\n', PreSalOutFile);
    Fig = figure;

    Ptitle = sprintf('Pre SAL: %s', CaseName);
    Xlabel = 'Radius (km)';
    Ylabel = 'Height (km)';

    % plot only the shading, not the contour lines
    contourf(R, Z, P_VT', 'LineStyle', 'none');
    set(gca, 'FontSize', Fsize);
    colorbar;

    hold on;
    [ Cup Hup ] = contour(R, Z, P_UP', UpClevs, 'LineColor', 'w');
    clabel(Cup, 'Color', 'w');

%    [ Cdn Hdn ] = contour(R, Z, P_DN', DnClevs, 'LineColor', 'w', 'LineStyle', '--');
%    clabel(Cdn, 'Color', 'w');

    title(Ptitle);
    xlabel(Xlabel);
    ylabel(Ylabel);

    caxis(VtClim);

    saveas(Fig, PreSalOutFile);
    close(Fig);


    % plot the SAL data
    fprintf('  Writing: %s\n', SalOutFile);
    Fig = figure;

    Ptitle = sprintf('SAL: %s', CaseName);
    Xlabel = 'Radius (km)';
    Ylabel = 'Height (km)';

    % plot only the shading, not the contour lines
    contourf(R, Z, S_VT', 'LineStyle', 'none');
    set(gca, 'FontSize', Fsize);
    colorbar;

    hold on;
    [ Cup Hup ] = contour(R, Z, S_UP', UpClevs, 'LineColor', 'w');
    clabel(Cup, 'Color', 'w');

%    [ Cdn Hdn ] = contour(R, Z, S_DN', DnClevs, 'LineColor', 'w', 'LineStyle', '--');
%    clabel(Cdn, 'Color', 'w');

    title(Ptitle);
    xlabel(Xlabel);
    ylabel(Ylabel);

    caxis(VtClim);

    saveas(Fig, SalOutFile);
    close(Fig);
    fprintf('\n');


  end
end
