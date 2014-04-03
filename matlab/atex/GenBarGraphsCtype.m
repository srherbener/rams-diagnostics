function [ ] = GenBarGraphsCtype(ConfigFile)
% GenBarGraphsCtype generate bar plots using cloud type filtered data

  % processing bgraph_<var>.h5 files
  % Two main variables: Averages and Npoints which are organized
  % as: (v,s,c,g) where
  %   v -> Variable name (cot_TSTART, eg)
  %   s -> SST   (ascending order)
  %   c -> CCN   (ascending order)
  %   g -> GCCN  (ascending order)

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);
    
  Ddir = Config.DiagDir;
  Pdir = Config.PlotDir;

  % Each entry in PlotDefs defines how to build a single plot. The format is:
  %  {
  %  Name
  %  Input file
  %  Title  --> { PanelMarker Label }
  %  Xlabel
  %  Ylabel
  %  LegSpec --> { LegText LegLoc }
  %  Bcolors
  %  Pstyle   ('grouped', 'stacked', etc.)
  %  VarSpec --> { Vname Select Transpose Xname Yname }
  %    Select has to be such that the result will reduce to a 2D array
  %       use vectors (okay to have 1 element vector) for selection
  %       eg: { [1] [1:3] [ 2 4 6 ] [3] } says to select
  %           Averages(1,1:3,[2 4 6],3) reducing it down to 2 dimensions
  %    Then Transpose says whether or not to transpose the array
  %    after it is selected down to two dimensions
  %       1 --> transpose
  %       0 --> do not transpose
  %  Output file
  %
  % Need to keep this in sync with GenBarGraphFilesCtype.m
  %   Var order (selection for the first dimension, v, of main variables)
  %     1  _TSTART
  %     2  _TMID
  %     3  _TEND
  %     4  _TALL
  %
  %     5  _strat_TSTART
  %     6  _strat_TMID
  %     7  _strat_TEND
  %     8  _strat_TALL
  %
  %     9  _scmix_TSTART
  %    10  _scmix_TMID
  %    11  _scmix_TEND
  %    12  _scmix_TALL
  %
  %    13  _scmix_TSTART
  %    14  _scmix_TMID
  %    15  _scmix_TEND
  %    16  _scmix_TALL

  PlotDefs = {
       % COT averages, all time points, CCN only
       {
       'COT Avg TALL CCN only'
       'bgraph_cot.h5'
       { 'a' 'All Points' }
       'N_a (# cm^-^3)'
       'COT'
       { [ 0 0 0 ] [ 0 1 1 ] [ 1 0 1 ] }
       { { '293 K', '298 K', '303 K' } 'NorthWest' }
       'grouped'
       { 'Averages' { [4] [1:3] [1:6] [1] } 1 'CCN' 'SST' }
       'bars_avg_cot_TALL_CO.jpg'
       }

       % PR averages, all time points, CCN only
       {
       'PR Avg TALL CCN only'
       'bgraph_pcprr.h5'
       { 'a' 'All Points' }
       'N_a (# cm^-^3)'
       'Precip Rate (mm h^-^1)'
       { [ 0 0 0 ] [ 0 1 1 ] [ 1 0 1 ] }
       { { '293 K', '298 K', '303 K' } 'NorthWest' }
       'grouped'
       { 'Averages' { [4] [1:3] [1:6] [1] } 1 'CCN' 'SST' }
       'bars_avg_pcprr_TALL_CO.jpg'
       }

       % LWP averages, all time points, CCN only
       {
       'LWP Avg TALL CCN only'
       'bgraph_lwp.h5'
       { 'a' 'All Points' }
       'N_a (# cm^-^3)'
       'LWP (mm)'
       { [ 0 0 0 ] [ 0 1 1 ] [ 1 0 1 ] }
       { { '293 K', '298 K', '303 K' } 'NorthWest' }
       'grouped'
       { 'Averages' { [4] [1:3] [1:6] [1] } 1 'CCN' 'SST' }
       'bars_avg_lwp_TALL_CO.jpg'
       }

       % Cloud depth averages, all time points, CCN only
       {
       'Cloud Depth Avg TALL CCN only'
       'bgraph_cdepth.h5'
       { 'a' 'All Points' }
       'N_a (# cm^-^3)'
       'Cloud Depth (m)'
       { [ 0 0 0 ] [ 0 1 1 ] [ 1 0 1 ] }
       { { '293 K', '298 K', '303 K' } 'NorthWest' }
       'grouped'
       { 'Averages' { [4] [1:3] [1:6] [1] } 1 'CCN' 'SST' }
       'bars_avg_cdepth_TALL_CO.jpg'
       }
     };


  for ipd = 1:length(PlotDefs)
    PdName  = PlotDefs{ipd}{1};
    InFile  = sprintf('%s/%s', Ddir, PlotDefs{ipd}{2});

    % Title
    Pmark   = PlotDefs{ipd}{3}{1};
    Ptitle  = PlotDefs{ipd}{3}{2};

    Xlabel  = PlotDefs{ipd}{4};
    Ylabel  = PlotDefs{ipd}{5};
    Bcolors = PlotDefs{ipd}{6};

    % LegSpec
    LegText = PlotDefs{ipd}{7}{1};
    LegLoc  = PlotDefs{ipd}{7}{2};

    Pstyle  = PlotDefs{ipd}{8};

    % VarSpec
    InVarName = PlotDefs{ipd}{9}{1};
    Vsel      = PlotDefs{ipd}{9}{2}{1};
    Ssel      = PlotDefs{ipd}{9}{2}{2};
    Csel      = PlotDefs{ipd}{9}{2}{3};
    Gsel      = PlotDefs{ipd}{9}{2}{4};
    Tpose     = PlotDefs{ipd}{9}{3};
    Xname     = PlotDefs{ipd}{9}{4};
    Yname     = PlotDefs{ipd}{9}{5};

    OutFname = PlotDefs{ipd}{10};

    fprintf('********************************************************************\n');
    fprintf('Generating bar graph: %s\n', PdName);
    fprintf('  Input File: %s\n', InFile);
    fprintf('    Variable: %s\n', InVarName);
    fprintf('  Selection specs:\n');
    fprintf('    Vsel: %s\n', strtrim(sprintf('%d ', Vsel)));
    fprintf('    Ssel: %s\n', strtrim(sprintf('%d ', Ssel)));
    fprintf('    Csel: %s\n', strtrim(sprintf('%d ', Csel)));
    fprintf('    Gsel: %s\n', strtrim(sprintf('%d ', Gsel)));
    fprintf('  Transpose: %d\n', Tpose);
    fprintf('  Xname: %s\n', Xname);
    fprintf('  Yname: %s\n', Yname);
    fprintf('\n');
 

    % read in variables
    HDATA = hdf5read(InFile, InVarName);
    BDATA = squeeze(HDATA(Vsel, Ssel, Csel, Gsel));
    if (Tpose == 1)
      BDATA = BDATA';
    end

    XVAR = hdf5read(InFile, Xname);
    YVAR = hdf5read(InFile, Yname);

    % do the plot
    clear AxisProps;
    iaxis = 0; % set this to next availble slot in the AxisProps array

    iaxis = iaxis + 1;
    AxisProps(iaxis).Name = 'FontSize';
    AxisProps(iaxis).Val  = 20;

    if (strcmp(Xname, 'CCN'))
      % label tick marks with CCN value
      Nvar = length(XVAR);
      X = 1:Nvar;      % use these for x-axis values -> evenly spaced integers creates evenly spaced bars
      clear XTIckLabels;
      for i = 1:Nvar
        XTickLabels{i} = sprintf('%d', XVAR(i));
      end
      iaxis = iaxis + 1;
      AxisProps(iaxis).Name = 'XTickLabel';
      AxisProps(iaxis).Val = XTickLabels;
    else
      X = XVAR;
    end

    % plot averages
    OutFile = sprintf('%s/%s', Pdir, OutFname);
    fprintf('      %s\n', OutFile);

    iaxis = iaxis + 1;
    AxisProps(iaxis).Name = 'Ylim';
    AxisProps(iaxis).Val  = [ 0 1.2*max(BDATA(:)) ];
    PlotBarSet(X, BDATA, Ptitle, Pmark, Xlabel, Ylabel, Bcolors, Pstyle, LegText, LegLoc, AxisProps, OutFile);

    fprintf('\n');

  end
end
