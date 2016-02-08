function [ ] = GenCintMeas()
% GenTsavgHists function to calculate temperature histograms

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
%    'TSD_SAL_NODUST'
%    'TSD_NONSAL_DUST'
%    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  % Description of measurements
  FileList = {
    { 'DIAGS/storm_profiles_<CASE>.h5' '/spath_i_aero_mass' 'Initial Non-Activated Dust' }
    { 'DIAGS/storm_profiles_<CASE>.h5' '/spath_f_aero_mass' 'Final Non-Activated Dust' }
    { 'DIAGS/storm_profiles_<CASE>.h5' '/spath_f_aero_mass_delta' 'Non-Activated Dust Difference' }

    { 'DIAGS/storm_profiles_<CASE>.h5' '/spath_f_dust_hydro' 'Final Activated Dust' }
    };
  Nfiles = length(FileList);

  Zlev = 5.5; % km

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating column integrated measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    for ifile = 1:Nfiles
      Ifile    = FileList{ifile}{1};
      Ivname   = FileList{ifile}{2};
      Quantity = FileList{ifile}{3};

      InFile  = regexprep(Ifile, '<CASE>', Case);

      fprintf('  Reading: %s (%s)\n', InFile, Ivname);

      % HDATA will be (x,y,t)
      HDATA = squeeze(h5read(InFile, Ivname));
      X = squeeze(h5read(InFile, '/x_coords'));
      Y = squeeze(h5read(InFile, '/y_coords'));
      Z = squeeze(h5read(InFile, '/z_coords'));
      T = squeeze(h5read(InFile, '/t_coords'));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      Z_KM = Z ./ 1000;
      Z1 = find(Z_KM >= Zlev, 1, 'first');

      % HDATA is (z), ie, a column vector
      % Integrate and report result
      % Just repeat delta z for top grid cell (should be the same given RAMS max delta z
      % configuration)
      DELTA_Z = Z(2:end) - Z(1:end-1);
      DELTA_Z(Nz) = DELTA_Z(Nz-1);
      ZINT = sum(HDATA .* DELTA_Z);

      ZINT_HLEV = sum(HDATA(Z1:end) .* DELTA_Z(Z1:end));

      % Convert from ug/m2 to ug/cm2
      ZINT = ZINT .* 1e-4;
      ZINT_HLEV = ZINT_HLEV .* 1e-4;
      
      fprintf('    Integrated %s: %.3f (ug cm^-^2)\n', Quantity, ZINT);
      fprintf('    Integrated upper level (Z > 6km) %s: %.3f (ug cm^-^2)\n', Quantity, ZINT_HLEV);
      fprintf('\n');
    end
  end
end 
