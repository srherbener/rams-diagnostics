function [ ] = GenVintTseries()
% GenVintTseries function to calculate time series of vertically integrated quantities

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
%    'TSD_SAL_NODUST'
%    'TSD_NONSAL_DUST'
%    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  Zlev = 6; % km

  % Description of measurements
  MeasList = {
    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_dust_hydro' '/spath_dust_hydro' }

    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_dust_cloud' '/spath_dust_cloud' }
    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_dust_rain'  '/spath_dust_rain'  }
    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_dust_pris'  '/spath_dust_pris'  }
    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_dust_snow'  '/spath_dust_snow'  }
    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_dust_aggr'  '/spath_dust_aggr'  }
    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_dust_graup' '/spath_dust_graup' }
    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_dust_hail'  '/spath_dust_hail'  }

    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_d1_mass'     '/spath_d1_mass'     }
    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_d2_mass'     '/spath_d2_mass'     }
    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_ccn_mass'    '/spath_ccn_mass'    }
    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_aero_mass'   '/spath_aero_mass'   }
    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_ra_mass'     '/spath_ra_mass'     }
    { 'DIAGS/storm_hovmollers_<CASE>.h5' '/spath_tracer_mass' '/spath_tracer_mass' }
    };
  Nmeas = length(MeasList);

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating vertically integrated measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    OutFile = sprintf('DIAGS/vint_meas_%s.h5', Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    for imeas = 1:Nmeas
      Ifile    = MeasList{imeas}{1};
      Ivname   = MeasList{imeas}{2};
      Ovname   = MeasList{imeas}{3};

      InFile  = regexprep(Ifile, '<CASE>', Case);
      fprintf('  Reading: %s (%s)\n', InFile, Ivname);

      % HDATA will be (z,t)
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

      % Integrate across z dimension
      % Just repeat delta z for top grid cell (should be the same given RAMS max delta z
      % configuration)
      DELTA_Z = Z(2:end) - Z(1:end-1);
      DELTA_Z(Nz) = DELTA_Z(Nz-1);
      DELTA_Z = repmat(DELTA_Z, [ 1 Nt ]);
      ZINT = squeeze(nansum((HDATA .* DELTA_Z), 1));

      ZINT_HLEV = squeeze(nansum((HDATA(Z1:end,:) .* DELTA_Z(Z1:end,:)), 1));

      % Convert from ug/m2 to ug/cm2
      ZINT      = ZINT .* 1e-4;
      ZINT_HLEV = ZINT_HLEV .* 1e-4;
      
      % If this is the first set, write out the coordinates into the output file
      % so that subsequent variables can have these coordinates attached.
      if (imeas == 1)
        Xname = '/x_coords';
        Yname = '/y_coords';
        Zname = '/z_coords';
        Tname = '/t_coords';
  
        CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
        % Add COARDS annotations
        NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
      end

      % ZINT will be a vector: (t)
      Vsize = Nt;
      DimOrder = { 't' };

      OvnameHlev = sprintf('%s_hlev', Ovname);

      fprintf('  Writing: %s (%s)\n', OutFile, Ovname);
      h5create(OutFile, Ovname, Vsize);
      h5write(OutFile, Ovname, ZINT);

      fprintf('  Writing: %s (%s)\n', OutFile, OvnameHlev);
      h5create(OutFile, OvnameHlev, Vsize);
      h5write(OutFile, OvnameHlev, ZINT_HLEV);

      AttachDimensionsXyzt(OutFile, Ovname,     DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, OvnameHlev, DimOrder, Xname, Yname, Zname, Tname);
    end
  end
end 
