function [ ] = GenDustAdvection()
% GenDustAdvection function to calculate advection at the edges of the analyses regions

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
%    'TSD_SAL_NODUST'
%    'TSD_NONSAL_DUST'
%    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  % Description of measurements
  %   Four 
  RegionList = {
    { 'FILTERS/sal_<CASE>.h5'      '/filter' 'sal'    } 

%    { 'FILTERS/sal_ar_<CASE>.h5'   '/filter' 'sal_ar' }
%    { 'FILTERS/sal_path_<CASE>.h5' '/filter' 'spath'  }
%    { 'FILTERS/all_500_<CASE>.h5'  '/filter' 'storm'  }
    };
  Nregions = length(RegionList);

  % data to calculate advection terms
  AdvList = {
     % momentum
    { 'HDF5/<CASE>/HDF5/u-<CASE>-AS-2006-08-20-120000-g3.h5' '/u' }
    { 'HDF5/<CASE>/HDF5/v-<CASE>-AS-2006-08-20-120000-g3.h5' '/v' }
    { 'HDF5/<CASE>/HDF5/z-<CASE>-AS-2006-08-20-120000-g3.h5' '/z' }

    % dust mass
    { 'HDF5/<CASE>/HDF5/dust_mass-<CASE>-AS-2006-08-20-120000-g3.h5'  '/dust_mass'       }
    { 'HDF5/<CASE>/HDF5/dust_hydro-<CASE>-AS-2006-08-20-120000-g3.h5' '/dust_hydro_mass' }
    { 'HDF5/<CASE>/HDF5/ra_mass-<CASE>-AS-2006-08-20-120000-g3.h5'    '/ra_mass'         }
    };

whos

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating dust advection measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Place all measurements into one case specific output file.
    % If the file exists, remove it so that the HDF5 commands
    % can effectively re-create datasets.
    OutFile = sprintf('DIAGS/residual_mass_%s.h5', Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    for ireg = 1:Nregions
      % The filter defines the edges of the horizontal domain where we want
      % to do the advection calculations.
      FilterFtemplate = RegionList{ireg}{1};
      FilterVname     = RegionList{ireg}{2};
      RegName         = RegionList{ireg}{3};

      FilterFile = regexprep(FilterFtemplate, '<CASE>', Case);
      fprintf('  Reading: %s (%s)\n', FilterFile, FilterVname);

      X = squeeze(h5read(FilterFile, '/x_coords'));
      Y = squeeze(h5read(FilterFile, '/y_coords'));
      Z = squeeze(h5read(FilterFile, '/z_coords'));
      T = squeeze(h5read(FilterFile, '/t_coords'));

      Nx       = length(X);
      Ny       = length(Y);
      FilterNz = length(Z);
      Nt       = length(T);

      % momentum and dust mass files
      Ufile = regexprep(AdvList{1}{1}, '<CASE>', Case);
      Uvname = AdvList{1}{2};
      fprintf('  Reading: %s (%s)\n', Ufile, Uvname);

      Vfile = regexprep(AdvList{2}{1}, '<CASE>', Case);
      Vvname = AdvList{2}{2};
      fprintf('  Reading: %s (%s)\n', Vfile, Vvname);

      Zfile = regexprep(AdvList{3}{1}, '<CASE>', Case);
      Zvname = AdvList{3}{2};
      fprintf('  Reading: %s (%s)\n', Zfile, Zvname);
      
      MdFile = regexprep(AdvList{4}{1}, '<CASE>', Case);
      MdVname = AdvList{4}{2};
      fprintf('  Reading: %s (%s)\n', MdFile, MdVname);

      MdhyFile = regexprep(AdvList{5}{1}, '<CASE>', Case);
      MdhyVname = AdvList{5}{2};
      fprintf('  Reading: %s (%s)\n', MdhyFile, MdhyVname);
      
      MdrgnFile = regexprep(AdvList{6}{1}, '<CASE>', Case);
      MdrgnVname = AdvList{6}{2};
      fprintf('  Reading: %s (%s)\n', MdrgnFile, MdrgnVname);

      Z = squeeze(h5read(Ufile, '/z_coords'));
      Nz = length(Z);

      fprintf('    Input variable dimensions:\n');
      fprintf('      Nx = %d\n', Nx);
      fprintf('      Ny = %d\n', Ny);
      fprintf('      Nz = %d\n', Nz);
      fprintf('      Nt = %d\n', Nt);
      fprintf('\n');
      fprintf('      FilterNz = %d\n', FilterNz);
      fprintf('\n');

      % process one time step at a time in order to keep memory usage down
      for it = 1:Nt
        FILTER = squeeze(h5read(FilterFile, FilterVname, [ 1 1 1 it ], [ Nx Ny FilterNz 1 ]));

        U = squeeze(h5read(Ufile, Uvname, [ 1 1 1 it ], [ Nx Ny Nz 1 ]));

        if (mod(it, 10) == 0)
          fprintf('  Completed time step %d out of %d\n', it, Nt);
        end
      end

  
%%%      % Total unactivated dust mass, all levels
%%%  
%%%      X = squeeze(h5read(InFile, '/x_coords'));
%%%      Y = squeeze(h5read(InFile, '/y_coords'));
%%%      Z = squeeze(h5read(InFile, '/z_coords'));
%%%      T = squeeze(h5read(InFile, '/t_coords'));
%%%  
%%%      Nx = length(X);
%%%      Ny = length(Y);
%%%      Nz = length(Z);
%%%      Nt = length(T);
%%%  
%%%      % regenerated dust mass, all levels
%%%      InVname = sprintf('/%s_ra_total_mass', Vprefix);
%%%      fprintf('  Reading: %s (%s)\n', InFile, InVname);
%%%      TS_MDRGN = squeeze(h5read(InFile, InVname));
%%%  
%%%      % dust in hydrometeor mass, all levels
%%%      InVname = sprintf('/%s_dust_hydro_total_mass', Vprefix);
%%%      fprintf('  Reading: %s (%s)\n', InFile, InVname);
%%%      TS_MDHY = squeeze(h5read(InFile, InVname));
%%%  
%%%      % dust deposited on surface
%%%      InVname = sprintf('/%s_dust_sfc_total_mass', Vprefix);
%%%      fprintf('  Reading: %s (%s)\n', InFile, InVname);
%%%      TS_MDSFC = squeeze(h5read(InFile, InVname));
%%%  
%%%      fprintf('\n');
%%%  
%%%  
%%%      % Total unactivated dust mass, upper levels
%%%      InVname = sprintf('/%s_dust_total_mass_hlev', Vprefix);
%%%      fprintf('  Reading: %s (%s)\n', InFile, InVname);
%%%      TS_MD_HLEV = squeeze(h5read(InFile, InVname));
%%%  
%%%      % regenerated dust mass, upper levels
%%%      InVname = sprintf('/%s_ra_total_mass_hlev', Vprefix);
%%%      fprintf('  Reading: %s (%s)\n', InFile, InVname);
%%%      TS_MDRGN_HLEV = squeeze(h5read(InFile, InVname));
%%%  
%%%      % dust in hydrometeor mass, upper levels
%%%      InVname = sprintf('/%s_dust_hydro_total_mass_hlev', Vprefix);
%%%      fprintf('  Reading: %s (%s)\n', InFile, InVname);
%%%      TS_MDHY_HLEV = squeeze(h5read(InFile, InVname));
%%%  
%%%      fprintf('\n');
%%%  
%%%      % residual is = Md - (Mdrgn + Mdhy + Mdsfc)
%%%      % Mdsfc is only for sfc to tropopause integrated values
%%%      TS_MDRES = TS_MD - (TS_MDRGN + TS_MDHY + TS_MDSFC);
%%%      TS_MDRES_HLEV = TS_MD_HLEV - (TS_MDRGN_HLEV + TS_MDHY_HLEV);
%%%
%%%      % MD_RMVD is time series of removed amount of dust (from original amount)
%%%      % residual removed is = MD_RMVD - (Mdrgn + Mdhy + Mdsfc)
%%%      % Mdsfc is only for sfc to tropopause integrated values
%%%      TS_MD_RMVD = TS_MD(1) - TS_MD;
%%%      TS_MD_RMVD_HLEV = TS_MD_HLEV(1) - TS_MD_HLEV;
%%%
%%%      TS_MDRES_RMVD = TS_MD_RMVD - (TS_MDRGN + TS_MDHY + TS_MDSFC);
%%%      TS_MDRES_RMVD_HLEV = TS_MD_RMVD_HLEV - (TS_MDRGN_HLEV + TS_MDHY_HLEV);
%%%
%%%      % Calculate rates: d(residual) / d(time)
%%%      TS_MDRES_RATE = (TS_MDRES(2:end) - TS_MDRES(1:end-1)) ./ (T(2:end) - T(1:end-1));
%%%      TS_MDRES_HLEV_RATE = (TS_MDRES_HLEV(2:end) - TS_MDRES_HLEV(1:end-1)) ./ (T(2:end) - T(1:end-1));
%%%
%%%      TS_MDRES_RMVD_RATE = (TS_MDRES_RMVD(2:end) - TS_MDRES_RMVD(1:end-1)) ./ (T(2:end) - T(1:end-1));
%%%      TS_MDRES_RMVD_HLEV_RATE = (TS_MDRES_RMVD_HLEV(2:end) - TS_MDRES_RMVD_HLEV(1:end-1)) ./ (T(2:end) - T(1:end-1));
%%%  
%%%      % Write out results
%%%      Vsize = Nt;
%%%      DimOrder = { 't' };
%%%  
%%%      if (ireg == 1)
%%%        % create coordinates
%%%        Xname = '/x_coords';
%%%        Yname = '/y_coords';
%%%        Zname = '/z_coords';
%%%        Tname = '/t_coords';
%%%  
%%%        CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
%%%        % Add COARDS annotations
%%%        NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
%%%      end
%%%  
%%%      % residual, all levels
%%%      OutVname = sprintf('/%s_residual_total_mass', Vprefix);
%%%      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
%%%      h5create(OutFile, OutVname, Vsize);
%%%      h5write(OutFile, OutVname, TS_MDRES);
%%%      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
%%%  
%%%      % residual, upper levels
%%%      OutVname = sprintf('/%s_residual_total_mass_hlev', Vprefix);
%%%      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
%%%      h5create(OutFile, OutVname, Vsize);
%%%      h5write(OutFile, OutVname, TS_MDRES_HLEV);
%%%      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
%%%  
%%%      % residual removed, all levels
%%%      OutVname = sprintf('/%s_residual_rmvd_total_mass', Vprefix);
%%%      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
%%%      h5create(OutFile, OutVname, Vsize);
%%%      h5write(OutFile, OutVname, TS_MDRES_RMVD);
%%%      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
%%%  
%%%      % residual, upper levels
%%%      OutVname = sprintf('/%s_residual_rmvd_total_mass_hlev', Vprefix);
%%%      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
%%%      h5create(OutFile, OutVname, Vsize);
%%%      h5write(OutFile, OutVname, TS_MDRES_RMVD_HLEV);
%%%      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
%%%
%%%      %%%%% RATES %%%%%%
%%%      Vsize = Nt-1;
%%%      % residual, all levels
%%%      OutVname = sprintf('/%s_residual_total_mass_rate', Vprefix);
%%%      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
%%%      h5create(OutFile, OutVname, Vsize);
%%%      h5write(OutFile, OutVname, TS_MDRES_RATE);
%%%      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
%%%  
%%%      % residual, upper levels
%%%      OutVname = sprintf('/%s_residual_total_mass_hlev_rate', Vprefix);
%%%      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
%%%      h5create(OutFile, OutVname, Vsize);
%%%      h5write(OutFile, OutVname, TS_MDRES_HLEV_RATE);
%%%      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
%%%  
%%%      % residual removed, all levels
%%%      OutVname = sprintf('/%s_residual_rmvd_total_mass_rate', Vprefix);
%%%      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
%%%      h5create(OutFile, OutVname, Vsize);
%%%      h5write(OutFile, OutVname, TS_MDRES_RMVD_RATE);
%%%      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
%%%  
%%%      % residual, upper levels
%%%      OutVname = sprintf('/%s_residual_rmvd_total_mass_hlev_rate', Vprefix);
%%%      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
%%%      h5create(OutFile, OutVname, Vsize);
%%%      h5write(OutFile, OutVname, TS_MDRES_RMVD_HLEV_RATE);
%%%      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
%%%
%%%      fprintf('\n');
    end
  end
end 
