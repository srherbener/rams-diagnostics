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
  InFileTemplate = 'DIAGS/storm_xsections_<CASE>.h5';
  PreSalVtVname  = '/ps_speed_t';
  SalVtVname     = '/s_speed_t';
  PreSalUpVname  = '/ps_updraft';
  SalUpVname     = '/s_updraft';
  PreSalDnVname  = '/ps_dndraft';
  SalDnVname     = '/s_dndraft';

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

    InFile = regexprep(InFileTemplate, '<CASE>', Case);

    PreSalOutFile = regexprep(PreSalOfileTemplate, '<CASE>', Case);
    SalOutFile    = regexprep(SalOfileTemplate, '<CASE>', Case);

    fprintf('  Reading: %s (%s, %s)\n', InFile, PreSalVtVname, SalVtVname);
    fprintf('  Reading: %s (%s, %s)\n', InFile, PreSalUpVname, SalUpVname);
    fprintf('  Reading: %s (%s, %s)\n', InFile, PreSalDnVname, SalDnVname);
    fprintf('\n');

    % Read in vars. Vars are organized as (r,z,t) where
    % r is the radius, z is the height, and t is time.
    P_VT = squeeze(h5read(InFile, PreSalVtVname));
    S_VT = squeeze(h5read(InFile, SalVtVname));
    P_UP = squeeze(h5read(InFile, PreSalUpVname));
    S_UP = squeeze(h5read(InFile, SalUpVname));
    P_DN = squeeze(h5read(InFile, PreSalDnVname));
    S_DN = squeeze(h5read(InFile, SalDnVname));
    R  = squeeze(h5read(InFile, '/x_coords')) ./ 1000; % radius in km
    Z  = squeeze(h5read(InFile, '/z_coords')) ./ 1000; % height in km
    T  = (squeeze(h5read(InFile, '/t_coords')) ./ 3600) - 42; % time in hours, 0h -> sim time 42h

    % Find indicies corresponding to selection
    R1 = find(R >= Rinside,  1, 'first');
    R2 = find(R <= Routside, 1, 'last');

    Z1 = find(Z >= Zbot, 1, 'first');
    Z2 = find(Z <= Ztop, 1, 'last');

    % Vars are (r,z)
    P_VT = P_VT(R1:R2,Z1:Z2);
    P_UP = P_UP(R1:R2,Z1:Z2);
    P_DN = P_DN(R1:R2,Z1:Z2);

    S_VT = S_VT(R1:R2,Z1:Z2);
    S_UP = S_UP(R1:R2,Z1:Z2);
    S_DN = S_DN(R1:R2,Z1:Z2);

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
