function [ ] = GenTotalMassMeas()
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
  % third argument is scaling factor to convert to g/m2 or g/m3
  FileList = {
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5'   '/sal_sum_dust_sfc'    1e-2 '/sal_dust_sfc_total_mass'    }  % units are ug/cm2 so scale by 1e-2 to get g/m2
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5'   '/sal_sum_dust_mass'   1e-6 '/sal_dust_total_mass'        }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5'   '/sal_sum_tracer_mass' 1e-6 '/sal_tracer_total_mass'      }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5'   '/sal_sum_dust_hydro'  1e-6 '/sal_dust_hydro_total_mass'  }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_aero_<CASE>.h5'   '/sal_sum_aero_mass'   1e-6 '/sal_aero_total_mass'        }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_ccn_<CASE>.h5'    '/sal_sum_ccn_mass'    1e-6 '/sal_ccn_total_mass'         }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_ra_<CASE>.h5'     '/sal_sum_ra_mass'     1e-6 '/sal_ra_total_mass'          }  % units are ug/m3 so scale by 1e-6 to get g/m3

    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5' '/spath_sum_dust_sfc'    1e-2 '/spath_dust_sfc_total_mass'    }  % units are ug/cm2 so scale by 1e-2 to get g/m2
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5' '/spath_sum_dust_mass'   1e-6 '/spath_dust_total_mass'        }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5' '/spath_sum_tracer_mass' 1e-6 '/spath_tracer_total_mass'      }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_dust_<CASE>.h5' '/spath_sum_dust_hydro'  1e-6 '/spath_dust_hydro_total_mass'  }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_aero_<CASE>.h5' '/spath_sum_aero_mass'   1e-6 '/spath_aero_total_mass'        }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_ccn_<CASE>.h5'  '/spath_sum_ccn_mass'    1e-6 '/spath_ccn_total_mass'         }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hda_meas_ts_ra_<CASE>.h5'   '/spath_sum_ra_mass'     1e-6 '/spath_ra_total_mass'          }  % units are ug/m3 so scale by 1e-6 to get g/m3

    { 'DIAGS/hist_sums_az_dust_<CASE>.h5' '/all_sum_dust_sfc'    1e-2 '/storm_dust_sfc_total_mass'    }  % units are ug/cm2 so scale by 1e-2 to get g/m2
    { 'DIAGS/hist_sums_az_dust_<CASE>.h5' '/all_sum_dust_mass'   1e-6 '/storm_dust_total_mass'        }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hist_sums_az_dust_<CASE>.h5' '/all_sum_tracer_mass' 1e-6 '/storm_tracer_total_mass'      }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hist_sums_az_dust_<CASE>.h5' '/all_sum_dust_hydro'  1e-6 '/storm_dust_hydro_total_mass'  }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hist_sums_az_aero_<CASE>.h5' '/all_sum_aero_mass'   1e-6 '/storm_aero_total_mass'        }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hist_sums_az_ccn_<CASE>.h5'  '/all_sum_ccn_mass'    1e-6 '/storm_ccn_total_mass'         }  % units are ug/m3 so scale by 1e-6 to get g/m3
    { 'DIAGS/hist_sums_az_ra_<CASE>.h5'   '/all_sum_ra_mass'     1e-6 '/storm_ra_total_mass'          }  % units are ug/m3 so scale by 1e-6 to get g/m3
    };
  Nfiles = length(FileList);

  Zlev = 6000; % m

  HorizArea = 9000 * 9000; % horizontal grid cell area, 3km X 3km

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Generating total mass measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Place all measurements into one case specific output file.
    % If the file exists, remove it so that the HDF5 commands
    % can effectively re-create datasets.
    OutFile = sprintf('DIAGS/total_mass_%s.h5', Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    for ifile = 1:Nfiles
      Ifile    = FileList{ifile}{1};
      Ivname   = FileList{ifile}{2};
      Scale    = FileList{ifile}{3};
      OutVname = FileList{ifile}{4};

      InFile  = regexprep(Ifile, '<CASE>', Case);

      fprintf('  Reading: %s (%s)\n', InFile, Ivname);

      % HDATA will be (x,y,t)
      HDATA = squeeze(h5read(InFile, Ivname)) .* Scale;
      X = squeeze(h5read(InFile, '/x_coords'));
      Y = squeeze(h5read(InFile, '/y_coords'));
      Z = squeeze(h5read(InFile, '/z_coords'));
      T = squeeze(h5read(InFile, '/t_coords'));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      Hsize = size(HDATA);

      % HDATA is either (z,t) or (t)
      Z1 = find(Z >= Zlev, 1, 'first');
      if ((Hsize(1) > 1) && (Hsize(2) > 1))
        % HDATA is (z,t)
        % Create a single column with delta z values
        %   Just repeat delta z for top grid cell (should be
        %   the same given RAMS max delta z configuration)
        DELTA_Z = Z(2:end) - Z(1:end-1);
        DELTA_Z(Nz) = DELTA_Z(Nz-1);

        % Repeat the column Nt times to create delta z matching
        % up with the size of HDATA.
        DELTA_Z = repmat(DELTA_Z, [ 1 Nt ]);

        % Do the vertical integration (summation)
        MEAS = sum(HDATA .* DELTA_Z);
        MEAS_HLEV = sum(HDATA(Z1:end,:) .* DELTA_Z(Z1:end,:));
      else
        MEAS = HDATA;
        MEAS_HLEV = HDATA; % dummy placeholder since no z structure
      end

      % Do the horizontal integration
      MEAS = MEAS .* HorizArea;
      MEAS_HLEV = MEAS_HLEV .* HorizArea;

      % If this is the first set, write out the coordinates into the output file
      % so that subsequent variables can have these coordinates attached.
      if (ifile == 1)
        Xname = '/x_coords';
        Yname = '/y_coords';
        Zname = '/z_coords';
        Tname = '/t_coords';
  
        CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
        % Add COARDS annotations
        NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
      end

      Vsize = Nt;
      DimOrder = { 't' };

      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
      h5create(OutFile, OutVname, Vsize);
      h5write(OutFile, OutVname, MEAS);

      OutHlevVname = sprintf('%s_hlev', OutVname);
      fprintf('  Writing: %s (%s)\n', OutFile, OutHlevVname);
      h5create(OutFile, OutHlevVname, Vsize);
      h5write(OutFile, OutHlevVname, MEAS_HLEV);

      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, OutHlevVname, DimOrder, Xname, Yname, Zname, Tname);

      fprintf('\n');
    end
  end
end 
