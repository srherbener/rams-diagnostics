function [ ] = GenResidualMassMeas()
% GenTotalMassMeas function to calculate total mass in volume or area

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
    { 'DIAGS/total_mass_<CASE>.h5' 'sal'    }
    { 'DIAGS/total_mass_<CASE>.h5' 'sal_ar' }
    { 'DIAGS/total_mass_<CASE>.h5' 'spath'  }
    { 'DIAGS/total_mass_<CASE>.h5' 'storm'  }
    };
  Nfiles = length(FileList);

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating residual mass measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Place all measurements into one case specific output file.
    % If the file exists, remove it so that the HDF5 commands
    % can effectively re-create datasets.
    OutFile = sprintf('DIAGS/residual_mass_%s.h5', Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    for ifile = 1:Nfiles
      Ftemplate = FileList{ifile}{1};
      Vprefix   = FileList{ifile}{2};

      InFile  = regexprep(Ftemplate, '<CASE>', Case);
  
      % Total unactivated dust mass, all levels
      InVname = sprintf('/%s_dust_total_mass', Vprefix);
      fprintf('  Reading: %s (%s)\n', InFile, InVname);
      TS_MD = squeeze(h5read(InFile, InVname));
  
      X = squeeze(h5read(InFile, '/x_coords'));
      Y = squeeze(h5read(InFile, '/y_coords'));
      Z = squeeze(h5read(InFile, '/z_coords'));
      T = squeeze(h5read(InFile, '/t_coords'));
  
      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);
  
      % regenerated dust mass, all levels
      InVname = sprintf('/%s_ra_total_mass', Vprefix);
      fprintf('  Reading: %s (%s)\n', InFile, InVname);
      TS_MDRGN = squeeze(h5read(InFile, InVname));
  
      % dust in hydrometeor mass, all levels
      InVname = sprintf('/%s_dust_hydro_total_mass', Vprefix);
      fprintf('  Reading: %s (%s)\n', InFile, InVname);
      TS_MDHY = squeeze(h5read(InFile, InVname));
  
      % dust deposited on surface
      InVname = sprintf('/%s_dust_sfc_total_mass', Vprefix);
      fprintf('  Reading: %s (%s)\n', InFile, InVname);
      TS_MDSFC = squeeze(h5read(InFile, InVname));
  
      fprintf('\n');
  
  
      % Total unactivated dust mass, upper levels
      InVname = sprintf('/%s_dust_total_mass_hlev', Vprefix);
      fprintf('  Reading: %s (%s)\n', InFile, InVname);
      TS_MD_HLEV = squeeze(h5read(InFile, InVname));
  
      % regenerated dust mass, upper levels
      InVname = sprintf('/%s_ra_total_mass_hlev', Vprefix);
      fprintf('  Reading: %s (%s)\n', InFile, InVname);
      TS_MDRGN_HLEV = squeeze(h5read(InFile, InVname));
  
      % dust in hydrometeor mass, upper levels
      InVname = sprintf('/%s_dust_hydro_total_mass_hlev', Vprefix);
      fprintf('  Reading: %s (%s)\n', InFile, InVname);
      TS_MDHY_HLEV = squeeze(h5read(InFile, InVname));
  
      fprintf('\n');
  
      % residual is = Md - (Mdrgn + Mdhy + Mdsfc)
      % Mdsfc is only for sfc to tropopause integrated values
      TS_MDRES = TS_MD - (TS_MDRGN + TS_MDHY + TS_MDSFC);
      TS_MDRES_HLEV = TS_MD_HLEV - (TS_MDRGN_HLEV + TS_MDHY_HLEV);
  
      % Write out results
      Vsize = Nt;
      DimOrder = { 't' };
  
      if (ifile == 1)
        % create coordinates
        Xname = '/x_coords';
        Yname = '/y_coords';
        Zname = '/z_coords';
        Tname = '/t_coords';
  
        CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
        % Add COARDS annotations
        NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
      end
  
      % residual, all levels
      OutVname = sprintf('/%s_residual_total_mass', Vprefix);
      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
      h5create(OutFile, OutVname, Vsize);
      h5write(OutFile, OutVname, TS_MDRES);
      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
  
      % residual, upper levels
      OutVname = sprintf('/%s_residual_total_mass_hlev', Vprefix);
      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
      h5create(OutFile, OutVname, Vsize);
      h5write(OutFile, OutVname, TS_MDRES_HLEV);
      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
  
      fprintf('\n');
    end
  end
end 
