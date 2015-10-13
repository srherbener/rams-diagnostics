function [ ] = GenPdfFilesCtype()
% GenPdfFilesCtype generate pdfs from cloud type hist data

  CaseList = {
    { 'z.atex.ccn0050.sst293' 12 36 }
%    { 'z.atex.ccn0100.sst293' 12 36 }
%    { 'z.atex.ccn0200.sst293' 12 36 }
%    { 'z.atex.ccn0400.sst293' 12 36 }
%    { 'z.atex.ccn0800.sst293' 12 36 }
    { 'z.atex.ccn1600.sst293' 12 36 }

    { 'z.atex.ccn0050.sst298' 12 36 }
%    { 'z.atex.ccn0100.sst298' 12 36 }
%    { 'z.atex.ccn0200.sst298' 12 36 }
%    { 'z.atex.ccn0400.sst298' 12 36 }
%    { 'z.atex.ccn0800.sst298' 12 36 }
    { 'z.atex.ccn1600.sst298' 12 36 }
    };
  Ncases = length(CaseList);

  OutFprefix = 'DIAGS/pdfs';

  InFprefix = 'hist_data_ctype';
  VarSets = {
    { 'Precip Rate' 'TsAveragedData/hist_pcprr_<CASE>.h5' '/hist_pcprr'  '/x_coords' '/pdf_pcprr'  }
    };
  Nvsets = length(VarSets);

  for icase = 1:Ncases
    Case   = CaseList{icase}{1};
    Tstart = CaseList{icase}{2};
    Tend   = CaseList{icase}{3};

    fprintf('***************************************************************\n');
    fprintf('Generating pdf data: %s\n', Case);
    fprintf('  Interval for time averaging: %f %f\n', Tstart, Tend);
    fprintf('\n');

    % Delete an existing file so that the subsequent writes won't collide with
    % existing datasets.
    OutFile = sprintf('%s_%s.h5', OutFprefix, Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    % Process all variables for this case.
    for ivset = 1:Nvsets
      Descrip     = VarSets{ivset}{1};
      InFtemplate = VarSets{ivset}{2};
      HistVname   = VarSets{ivset}{3};
      BinsVname   = VarSets{ivset}{4};
      PdfVname    = VarSets{ivset}{5};

      TavgPdfVname = sprintf('%s_tavg', PdfVname);
      TavgHistVname = sprintf('%s_tavg', HistVname);

      fprintf('  Variable: %s\n', Descrip);

      % Hist data will be (b,t) where b is number of bins.
      InFile = regexprep(InFtemplate, '<CASE>', Case);
      fprintf('    Reading: %s (%s)\n', InFile, HistVname);
      HIST = squeeze(h5read(InFile, HistVname));
      B    = squeeze(h5read(InFile, BinsVname));
      T    = squeeze(h5read(InFile, '/t_coords')) ./ 3600; % hours

      Nb = length(B);
      Nt = length(T);

      % Select the data for time averaging
      T1 = find(T >= Tstart, 1, 'first');
      T2 = find(T <= Tend,   1, 'last');
      HIST_TAVG = squeeze(sum(HIST(:,T1:T2), 2));

      % Form pdfs, for the time series (b,t) pdf, need to
      % replicate the sum of each histogram along the bins
      % dimension so that an elementwise divide can be performed.
      SUM_HIST = squeeze(sum(HIST, 1));
      SUM_HIST = repmat(SUM_HIST, [ Nb 1 ]);
      PDF = HIST ./ SUM_HIST;

      % Time averaged pdf. sum(HIST_TAVG) reduces to a single value
      % which will work as is in the elementwise divide
      PDF_TAVG = HIST_TAVG ./ sum(HIST_TAVG);

      % Write both the histograms and pdfs into the output file.
      fprintf('    Writing: %s (%s)\n', OutFile, HistVname);
      h5create(OutFile, HistVname, size(HIST));
      h5write (OutFile, HistVname, HIST);

      fprintf('    Writing: %s (%s)\n', OutFile, TavgHistVname);
      h5create(OutFile, TavgHistVname, size(HIST_TAVG));
      h5write (OutFile, TavgHistVname, HIST_TAVG);

      fprintf('    Writing: %s (%s)\n', OutFile, PdfVname);
      h5create(OutFile, PdfVname, size(PDF));
      h5write (OutFile, PdfVname, PDF);

      fprintf('    Writing: %s (%s)\n', OutFile, TavgPdfVname);
      h5create(OutFile, TavgPdfVname, size(PDF_TAVG));
      h5write (OutFile, TavgPdfVname, PDF_TAVG);

      % Write out the coordinate values
      h5create(OutFile, '/Bins', size(B));
      h5write (OutFile, '/Bins', B);

      h5create(OutFile, '/Time', size(T));
      h5write (OutFile, '/Time', T);
      
      fprintf('\n');
    end
  end
end
