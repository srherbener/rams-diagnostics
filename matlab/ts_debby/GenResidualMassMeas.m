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

    InFile  = sprintf('DIAGS/total_mass_%s.h5', Case);

    % Total unactivated dust mass, all levels
    InVname = '/sal_dust_total_mass';
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
    InVname = '/sal_ra_total_mass';
    fprintf('  Reading: %s (%s)\n', InFile, InVname);
    TS_MDRGN = squeeze(h5read(InFile, InVname));

    % dust in hydrometeor mass, all levels
    InVname = '/sal_dust_hydro_total_mass';
    fprintf('  Reading: %s (%s)\n', InFile, InVname);
    TS_MDHY = squeeze(h5read(InFile, InVname));

    % dust deposited on surface
    InVname = '/sal_dust_sfc_total_mass';
    fprintf('  Reading: %s (%s)\n', InFile, InVname);
    TS_MDSFC = squeeze(h5read(InFile, InVname));

    fprintf('\n');


    % Total unactivated dust mass, upper levels
    InVname = '/sal_dust_total_mass_hlev';
    fprintf('  Reading: %s (%s)\n', InFile, InVname);
    TS_MD_HLEV = squeeze(h5read(InFile, InVname));

    % regenerated dust mass, upper levels
    InVname = '/sal_ra_total_mass_hlev';
    fprintf('  Reading: %s (%s)\n', InFile, InVname);
    TS_MDRGN_HLEV = squeeze(h5read(InFile, InVname));

    % dust in hydrometeor mass, upper levels
    InVname = '/sal_dust_hydro_total_mass_hlev';
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

    % create coordinates
    Xname = '/x_coords';
    Yname = '/y_coords';
    Zname = '/z_coords';
    Tname = '/t_coords';

    CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
    % Add COARDS annotations
    NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);

    % residual, all levels
    OutVname = '/sal_residual_total_mass';
    fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
    h5create(OutFile, OutVname, Vsize);
    h5write(OutFile, OutVname, TS_MDRES);
    AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);

    % residual, upper levels
    OutVname = '/sal_residual_total_mass_hlev';
    fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
    h5create(OutFile, OutVname, Vsize);
    h5write(OutFile, OutVname, TS_MDRES_HLEV);
    AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);

    fprintf('\n');
  end
end 
