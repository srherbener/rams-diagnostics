function [ ] = ReportDpMeas()
% ReportDpMeas function to report measurements for the DP paper

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
    'TSD_SD'
    'TSD_SD_1G'
%    'TSD_SAL_NODUST'
%    'TSD_NONSAL_DUST'
%    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  % Description of measurements
  MeasList = {
%    { 'DIAGS/total_mass_<CASE>.h5' '/sal'    'SAL_AR_L' }
    { 'DIAGS/total_mass_<CASE>.h5' '/sal_ar' 'SAL_AR_M' }
%    { 'DIAGS/total_mass_<CASE>.h5' '/spath'  'SPATH'    }

    };
  Nmeas = length(MeasList);

  for icase = 1:Ncases
    Case = CaseList{icase};

    fprintf('*****************************************************************\n');
    fprintf('DP Measurements:\n');
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Run through entire file list first, record the numbers and dump out
    % the report last. This way the report is all in one section instead of
    % interspersed between the "Reading: ..." messages.

    clear RDATA;
    for imeas = 1:Nmeas
      IfilePfx = MeasList{imeas}{1};
      IvarPfx  = MeasList{imeas}{2};
      Ilabel   = MeasList{imeas}{3};

      % Data will be (t), and units are grams
      InFile  = regexprep(IfilePfx, '<CASE>', Case);
      
      % Md all levels
      InVname = sprintf('%s_dust_total_mass', IvarPfx);
      fprintf('    Reading: %s (%s)\n', InFile, InVname);
      MD = squeeze(h5read(InFile, InVname));

      % Mdsfc all levels
      InVname = sprintf('%s_dust_sfc_total_mass', IvarPfx);
      fprintf('    Reading: %s (%s)\n', InFile, InVname);
      MDSFC = squeeze(h5read(InFile, InVname));

      % Mdhy upper levels
      InVname = sprintf('%s_dust_hydro_total_mass_hlev', IvarPfx);
      fprintf('    Reading: %s (%s)\n', InFile, InVname);
      MDHY = squeeze(h5read(InFile, InVname));

      % Mdrgn upper levels
      InVname = sprintf('%s_ra_total_mass_hlev', IvarPfx);
      fprintf('    Reading: %s (%s)\n', InFile, InVname);
      MDRGN = squeeze(h5read(InFile, InVname));

      %     Quantity             How to calculate
      %
      %   Initial dust                MD(1)
      %
      %   Dust removed                MDSFC(end)
      %
      %   Dust transported            Find index of max(MDRGN) -> MaxInd
      %      to upper levels          Formulas: 
      %                                  Want 2nd number to be value at t = 30 min.
      %                                  This is array element 2 for TSD_SAL_DUST
      %                                  and linear interp. of first two elements
      %                                  for TSD_SD and TSD_SD_1G.
      %                                  MDHY(MaxInd)  - MDHY(t=30m)
      %                                  MDRGN(MaxInd) - MDRGN(t=30m)
      %
      % Use max of MDRGN since this is the dust that will linger in the
      % upper levels after the storm passes. Also because MDRGN tends to
      % be the greater of the two. Convert to Tg.
      Minit = MD(1) * 1e-12;
      Mrmvd = MDSFC(end) * 1e-12;
     
      switch(Case)
        case { 'TSD_SAL_DUST' }
          Mdrgn1 = MDRGN(2);
          Mdhy1  = MDHY(2);

        case { 'TSD_SD' 'TSD_SD_1G' }
          Mdrgn1 = (MDRGN(1) + MDRGN(2)) * 0.5;
          Mdhy1  = (MDHY(1) + MDHY(2)) * 0.5;
      end
      Mdrgn1 = Mdrgn1 * 1e-12;
      Mdhy1  = Mdhy1 * 1e-12;

      [ Mdrgn2 MaxInd ] = max(MDRGN);
      Mdrgn2 = Mdrgn2 * 1e-12;
      Mdhy2  = MDHY(MaxInd) * 1e-12;

      % record for report
      RDATA{imeas} = { Ilabel Minit Mrmvd Mdhy1 Mdhy2 Mdrgn1 Mdrgn2 };

      fprintf('\n');
    end

    % Dump out report
    OutFile = sprintf('DIAGS/DpMeasReport_%s.txt', Case);
    fprintf('    Writing: %s\n', OutFile);

    OutFid = fopen(OutFile, 'w');

    fprintf(OutFid, '%-10s %10s %10s %10s %10s %10s %10s %10s %10s\n', 'Analysis', 'Md',      'Mdsfc', '', 'Mdhy',   '', 'Mdrgn',  '', 'Total');
    fprintf(OutFid, '%-10s %10s %10s %10s %10s %10s %10s %10s %10s\n', 'Region',   'Initial', 'Final', '', 'Lofted', '', 'Lofted', '', 'Lofted');
    fprintf(OutFid, '%-10s %10s %10s %10s %10s %10s %10s %10s %10s\n', '',         '(Tg)',    '(Tg)',  '', '(Tg)',   '', '(Tg)',   '', '(Tg)');
    fprintf(OutFid, '\n');
    for imeas = 1:length(RDATA)
      Ilabel     = RDATA{imeas}{1};
      Minit      = RDATA{imeas}{2};
      Mrmvd      = RDATA{imeas}{3};
      Mdhy1      = RDATA{imeas}{4};
      Mdhy2      = RDATA{imeas}{5};
      Mdrgn1     = RDATA{imeas}{6};
      Mdrgn2     = RDATA{imeas}{7};

      MdrgnLoft = Mdrgn2 - Mdrgn1;
      MdhyLoft  = Mdhy2  - Mdhy1;
      Mloft     = MdhyLoft + MdrgnLoft;

      MrmvdPct = sprintf('(%7.4f%%)', 100 * (Mrmvd/Minit));
      MdhyPct  = sprintf('(%7.4f%%)', 100 * (MdhyLoft/Minit));
      MdrgnPct = sprintf('(%7.4f%%)', 100 * (MdrgnLoft/Minit));
      MloftPct = sprintf('(%7.4f%%)', 100 * (Mloft/Minit));
      
      fprintf(OutFid, '%-10s %10.2e %10.2e %10s %10.2e %10s %10.2e %10s %10.2e %10s\n', ...
        Ilabel, Minit, Mrmvd, MrmvdPct, MdhyLoft, MdhyPct, MdrgnLoft, MdrgnPct, Mloft, MloftPct);

    end
    fclose(OutFid);
    fprintf('\n');
  end
end 
