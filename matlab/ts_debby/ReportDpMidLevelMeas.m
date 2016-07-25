function [ ] = ReportDpMidLevelMeas()
% ReportDpMeas function to report measurements for the DP paper

  % Cases
  CaseList = {
    'TSD_SAL_DUST'
%    'TSD_SAL_NODUST'
%    'TSD_NONSAL_DUST'
%    'TSD_NONSAL_NODUST'
    };
  Ncases = length(CaseList);

  % Description of measurements
  MeasList = {
    { 'DIAGS/total_mass_<CASE>.h5' '/sal_ar' 'SAL_AR_M' }

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

      % Mdrgn mid levels
      InVname = sprintf('%s_ra_total_mass_mlev', IvarPfx);
      fprintf('    Reading: %s (%s)\n', InFile, InVname);
      MDRGN = squeeze(h5read(InFile, InVname));

      % Mdrgn mid levels
      InVname = sprintf('%s_ra1_total_mass_mlev', IvarPfx);
      fprintf('    Reading: %s (%s)\n', InFile, InVname);
      MDRGN1 = squeeze(h5read(InFile, InVname));

      % Mdrgn mid levels
      InVname = sprintf('%s_ra2_total_mass_mlev', IvarPfx);
      fprintf('    Reading: %s (%s)\n', InFile, InVname);
      MDRGN2 = squeeze(h5read(InFile, InVname));

      %     Quantity             How to calculate
      %
      %   Initial dust                MD(1)
      %
      % Use max of MDRGN since this is the dust that will linger in the
      % upper levels after the storm passes. Also because MDRGN tends to
      % be the greater of the two. Convert to Tg.
      Minit = MD(1) * 1e-12;
     
      MdrgnStart  = MDRGN(2) * 1e-12;
      Mdrgn1Start = MDRGN1(2) * 1e-12;
      Mdrgn2Start = MDRGN2(2) * 1e-12;

      [ MdrgnEnd MaxInd ] = max(MDRGN);
      MdrgnEnd = MdrgnEnd * 1e-12;
      Mdrgn1End = MDRGN1(MaxInd) * 1e-12;
      Mdrgn2End = MDRGN2(MaxInd) * 1e-12;

      % record for report
      RDATA{imeas} = { Ilabel Minit MdrgnStart Mdrgn1Start Mdrgn2Start MdrgnEnd Mdrgn1End Mdrgn2End };

      fprintf('\n');
    end

    % Dump out report
    OutFile = sprintf('DIAGS/DpMidLevelMeasReport_%s.txt', Case);
    fprintf('    Writing: %s\n', OutFile);

    OutFid = fopen(OutFile, 'w');

    fprintf(OutFid, '%-10s %10s %10s %10s %10s %10s %10s\n', 'Analysis', 'Md',      'Mdrgn', '', 'Mdrgn1',   '', 'Mdrgn2');
    fprintf(OutFid, '%-10s %10s %10s %10s %10s %10s %10s\n', 'Region',   'Initial', 'MidLev', '', 'MidLev', '', 'MidLev');
    fprintf(OutFid, '%-10s %10s %10s %10s %10s %10s %10s\n', '',         '(Tg)',    '(Tg)',  '', '(Tg)',   '', '(Tg)');
    fprintf(OutFid, '\n');
    for imeas = 1:length(RDATA)
      Ilabel      = RDATA{imeas}{1};
      Minit       = RDATA{imeas}{2};
      MdrgnStart  = RDATA{imeas}{3};
      Mdrgn1Start = RDATA{imeas}{4};
      Mdrgn2Start = RDATA{imeas}{5};
      MdrgnEnd    = RDATA{imeas}{6};
      Mdrgn1End   = RDATA{imeas}{7};
      Mdrgn2End   = RDATA{imeas}{8};

      MdrgnMidLev = MdrgnEnd - MdrgnStart;
      Mdrgn1MidLev = Mdrgn1End - Mdrgn1Start;
      Mdrgn2MidLev = Mdrgn2End - Mdrgn2Start;

      MdrgnPct  = sprintf('(%7.4f%%)', 100 * (MdrgnMidLev/Minit));
      Mdrgn1Pct = sprintf('(%7.4f%%)', 100 * (Mdrgn1MidLev/Minit));
      Mdrgn2Pct = sprintf('(%7.4f%%)', 100 * (Mdrgn2MidLev/Minit));
      
      fprintf(OutFid, '%-10s %10.2e %10.2e %10s %10.2e %10s %10.2e %10s\n', ...
        Ilabel, Minit, MdrgnMidLev, MdrgnPct, Mdrgn1MidLev, Mdrgn1Pct, Mdrgn2MidLev, Mdrgn2Pct);

    end
    fclose(OutFid);
    fprintf('\n');
  end
end 
