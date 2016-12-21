function [ ] = ReportFsMeas()
% ReportFsMeas function to report measurements for the FS paper

  % Description of measurements
  % <blank> says to insert a blank line in the report
  MeasList = {
    'avg_wind'
    'ps_avg_wind'
    's_avg_wind'
    '<blank>'
    'max_wind'
    'ps_max_wind'
    's_max_wind'
    '<blank>'
    'avg_wind_t'
    'ps_avg_wind_t'
    's_avg_wind_t'
    '<blank>'
    'max_wind_t'
    'ps_max_wind_t'
    's_max_wind_t'
    };
  Nmeas = length(MeasList);

  InFile = 'DIAGS/fs_factors.h5';

  fprintf('*****************************************************************\n');
  fprintf('DP Measurements:\n');
  fprintf('\n');

  % Run through entire file list first, record the numbers and dump out
  % the report last. This way the report is all in one section instead of
  % interspersed between the "Reading: ..." messages.
  for imeas = 1:Nmeas
    Mname = MeasList{imeas};

    if (strcmp(Mname, '<blank>'))
      RDATA{imeas} = { Mname 0 0 0 0 0 0 0 };
    else
      AvgVname = sprintf('/%s_averages', Mname);
      FacVname = sprintf('/%s_factors', Mname);
  
      fprintf('  Reading: %s (%s, %s)\n', InFile, AvgVname, FacVname);
      AVGS = squeeze(h5read(InFile, AvgVname));
      FACS = squeeze(h5read(InFile, FacVname));

      % Data: name NSND SND NSD SD F1 F2 F12
      RDATA{imeas} = { Mname AVGS(1) AVGS(2) AVGS(3) AVGS(4) FACS(2) FACS(3) FACS(4) };
    end
  end
  fprintf('\n');

  % Dump out report
  OutFile = 'DIAGS/FsMeasReport.txt';
  fprintf('  Writing: %s\n', OutFile);

  OutFid = fopen(OutFile, 'w');

  fprintf(OutFid, '%15s     %10s %10s %10s %10s %10s %10s %10s\n', ...
    'Quantity', 'NSND',  'SND',   'NSD',   'SD',    'F1',    'F2',    'F12');
  fprintf(OutFid, '%15s     %10s %10s %10s %10s %10s %10s %10s\n', ...
    '',         '(m/s)', '(m/s)', '(m/s)', '(m/s)', '(m/s)', '(m/s)', '(m/s)');
  fprintf(OutFid, '\n');
  for imeas = 1:length(RDATA)
    Mname = RDATA{imeas}{1};
    NSND  = RDATA{imeas}{2};
    SND   = RDATA{imeas}{3};
    NSD   = RDATA{imeas}{4};
    SD    = RDATA{imeas}{5};
    F1    = RDATA{imeas}{6};
    F2    = RDATA{imeas}{7};
    F12   = RDATA{imeas}{8};

    if (strcmp(Mname, '<blank>'))
      fprintf(OutFid, '\n');
    else
      fprintf(OutFid, '%15s     %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', ...
        Mname, NSND, SND, NSD, SD, F1, F2, F12);
    end
  end
  fclose(OutFid);
  fprintf('\n');
end 
