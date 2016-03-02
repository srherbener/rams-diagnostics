function [ ] = ReportIntMassNumbers()
% ReportIntMassNumbers function to report initial and final volume integrated mass numbers

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
    { 'DIAGS/total_mass_<CASE>.h5'    '/sal_dust_total_mass'       'SAL_AR (large), Unactivated:'          }
    { 'DIAGS/total_mass_<CASE>.h5'    '/sal_dust_sfc_total_mass'   'SAL_AR (large), Deposited to surface:' }
    { 'DIAGS/total_mass_<CASE>.h5'    '/sal_dust_hydro_total_mass' 'SAL_AR (large), Inside hydrometeors:'  }
    { 'DIAGS/total_mass_<CASE>.h5'    '/sal_ra_total_mass'         'SAL_AR (large), Regenerated:'          }
    { 'DIAGS/residual_mass_<CASE>.h5' '/sal_residual_total_mass'   'SAL_AR (large), Advected (residual):'  }

    { 'DIAGS/total_mass_<CASE>.h5'    '/sal_ar_dust_total_mass'       'SAL_AR (med), Unactivated:'          }
    { 'DIAGS/total_mass_<CASE>.h5'    '/sal_ar_dust_sfc_total_mass'   'SAL_AR (med), Deposited to surface:' }
    { 'DIAGS/total_mass_<CASE>.h5'    '/sal_ar_dust_hydro_total_mass' 'SAL_AR (med), Inside hydrometeors:'  }
    { 'DIAGS/total_mass_<CASE>.h5'    '/sal_ar_ra_total_mass'         'SAL_AR (med), Regenerated:'          }
    { 'DIAGS/residual_mass_<CASE>.h5' '/sal_ar_residual_total_mass'   'SAL_AR (med), Advected (residual):'  }

    { 'DIAGS/total_mass_<CASE>.h5'    '/spath_dust_total_mass'       'SPATH (small), Unactivated:'          }
    { 'DIAGS/total_mass_<CASE>.h5'    '/spath_dust_sfc_total_mass'   'SPATH (small), Deposited to surface:' }
    { 'DIAGS/total_mass_<CASE>.h5'    '/spath_dust_hydro_total_mass' 'SPATH (small), Inside hydrometeors:'  }
    { 'DIAGS/total_mass_<CASE>.h5'    '/spath_ra_total_mass'         'SPATH (small), Regenerated:'          }
    { 'DIAGS/residual_mass_<CASE>.h5' '/spath_residual_total_mass'   'SPATH (small), Advected (residual):'  }

    { 'DIAGS/total_mass_<CASE>.h5'    '/storm_dust_total_mass'       'STORM (500km), Unactivated:'          }
    { 'DIAGS/total_mass_<CASE>.h5'    '/storm_dust_sfc_total_mass'   'STORM (500km), Deposited to surface:' }
    { 'DIAGS/total_mass_<CASE>.h5'    '/storm_dust_hydro_total_mass' 'STORM (500km), Inside hydrometeors:'  }
    { 'DIAGS/total_mass_<CASE>.h5'    '/storm_ra_total_mass'         'STORM (500km), Regenerated:'          }
    { 'DIAGS/residual_mass_<CASE>.h5' '/storm_residual_total_mass'   'STORM (500km), Advected (residual):'  }
    };
  Nfiles = length(FileList);
  LineSpace = 5; % put in a line space after every fifth line - matches grouping in FileList

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('Integrated Mass Measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Run through entire file list first, record the numbers and dump out
    % the report last. This way the report is all in one section instead of
    % interspersed between the "Reading: ..." messages.

    clear RDATA;
    for ifile = 1:Nfiles
      Ifile   = FileList{ifile}{1};
      InVname = FileList{ifile}{2};
      Descrip = FileList{ifile}{3};

      % HDATA will be (t), and units are grams
      InFile  = regexprep(Ifile, '<CASE>', Case);
      fprintf('  Reading: %s (%s)\n', InFile, InVname);
      HDATA = squeeze(h5read(InFile, InVname));

      Mstart = HDATA(1)   * 1e-12; % Convert to Tg
      Mend   = HDATA(end) * 1e-12; % Convert to Tg

      Mdiff = Mend - Mstart;

      % record the description, start value, end value and difference between start and end
      RDATA{ifile} = { Descrip Mstart Mend Mdiff };
    end
    fprintf('\n');

    % Dump out report
    OutFile = sprintf('DIAGS/IntMassReport_%s.txt', Case);
    fprintf('  Writing: %s\n', OutFile);

    OutFid = fopen(OutFile, 'w');

    fprintf(OutFid, '%-40s %10s %10s %10s\n', '        Mass Quantity', 'Start', 'End', 'Diff');
    fprintf(OutFid, '\n');
    for imeas = 1:length(RDATA)
      Descrip = RDATA{imeas}{1};
      Mstart  = RDATA{imeas}{2};
      Mend    = RDATA{imeas}{3};
      Mdiff   = RDATA{imeas}{4};
      fprintf(OutFid, '  %-40s %10f %10f %10f \n', Descrip, Mstart, Mend, Mdiff);

      if (mod(imeas, LineSpace) == 0)
        fprintf(OutFid, '\n');
      end
    end
    fclose(OutFid);
    fprintf('\n');
  end
end 
