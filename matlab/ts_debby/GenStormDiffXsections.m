function [ ] = GenStormDiffXsections()

  Ddir = 'DIAGS';
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % Description of diff pairs
  DiffList = {
    { '^all_ps' 'TSD_NONSAL_DUST'   'NSD'  'TSD_SAL_DUST'   'SD'  }
    { '^all_ps' 'TSD_NONSAL_NODUST' 'NSND' 'TSD_SAL_NODUST' 'SND' }

    { '^all_s' 'TSD_SAL_NODUST'    'SND'  'TSD_SAL_DUST'    'SD'  }
    { '^all_s' 'TSD_NONSAL_NODUST' 'NSND' 'TSD_NONSAL_DUST' 'NSD' }
    };

  Ndiffs = length(DiffList);

  InFileTemplate = 'DIAGS/storm_xsections_<CASE>.h5';
  OutFileTemplate = 'DIAGS/storm_diff_xsections_<CASE1>_<CASE2>.h5';

  for idiff = 1:Ndiffs
    Vselect    = DiffList{idiff}{1};
    Case1      = DiffList{idiff}{2};
    Case1Label = DiffList{idiff}{3};
    Case2      = DiffList{idiff}{4};
    Case2Label = DiffList{idiff}{5};

    fprintf('**********************************************************\n');
    fprintf('Generating storm difference cross-sections for cases: %s %s\n', Case1Label, Case2Label);
    fprintf('\n');

    InFile1 = regexprep(InFileTemplate, '<CASE>', Case1);
    InFile2 = regexprep(InFileTemplate, '<CASE>', Case2);

    fprintf('  Reading: %s\n', InFile1);
    fprintf('  Reading: %s\n', InFile2);
    fprintf('\n');

    % Place all vars into one case specific output file.
    % If the file exists, remove it so that the HDF5 commands
    % can effectively re-create datasets.
    OutFile = regexprep(OutFileTemplate, '<CASE1>', Case1Label);
    OutFile = regexprep(OutFile,         '<CASE2>', Case2Label);

    fprintf('  Writing: %s\n', InFile2);

    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    % Look in InFile1 to get a list of variables, assume that both input files have the
    % same list of variables.
    InFileInfo = h5info(InFile1);
    VarList = InFileInfo.Datasets;
    Nvars = length(VarList);

    % Get the coordinate values. Assume both input files have the same coordinate
    % descriptions.
    X = squeeze(h5read(InFile1, '/x_coords'));
    Y = squeeze(h5read(InFile1, '/y_coords'));
    Z = squeeze(h5read(InFile1, '/z_coords'));
    T = squeeze(h5read(InFile1, '/t_coords'));

    Nx = length(X);
    Ny = length(Y);
    Nz = length(Z);
    Nt = length(T);

    % Write out the coordinates into the output file so that subsequent variables
    % can have these coordinates attached.
    Xname = '/x_coords';
    Yname = '/y_coords';
    Zname = '/z_coords';
    Tname = '/t_coords';

    CreateDimensionsXyzt(OutFile, X, Y, Z, [ 1 ], Xname, Yname, Zname, Tname);
    % Add COARDS annotations
    NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);

    % Grab the datasets that match the patterns given Vselect
    for ivar = 1:Nvars
      Vname = VarList(ivar).Name;

      % skip if name doesn't match Vselect
      if (isempty(regexp(Vname, Vselect)))
        continue
      end

      % show only what we are processing
      fprintf('    Var: %s\n', Vname);

      InVname = sprintf('/%s', Vname);
      VAR1 = squeeze(h5read(InFile1, InVname));
      VAR2 = squeeze(h5read(InFile2, InVname));

      DIFF = VAR1 - VAR2;

      % Determine size of output var
      Vsize = size(DIFF);
      if ((Vsize(1) == 1) && (Vsize(2) == 1))
        % VAR is [ 1 1 ], ie. a scalar value
        Vsize = 1;
        DimOrder = { };
      elseif ((Vsize(1) == 1) || (Vsize(2) == 1))
        % VAR is [ 1 n ] or [ n 1 ], ie. a vector value
        % This implies 1D in radius which is the x dimension
        Vsize = Nx;
        DimOrder = { 'x' };
      else
        % Var is [ n m ]
        % Vsize does not need to be modified
        % This implies 2D in radius and height which are the x and z dimensions
        DimOrder = { 'x' 'z' };
      end
      OutVname = sprintf('%s_diff', InVname);

      % Collect output vars into one case specific file.
      h5create(OutFile, OutVname, Vsize);
      h5write(OutFile, OutVname, DIFF);

      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
    end
    fprintf('\n');

  end
end
