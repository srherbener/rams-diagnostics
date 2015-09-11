function [ ] = GenStormXsections()

  Ddir = 'DIAGS';
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  Cases = {
   'TSD_SAL_DUST'
   'TSD_SAL_NODUST'
   'TSD_NONSAL_DUST'
   'TSD_NONSAL_NODUST'
   };
  Ncases = length(Cases);

  % Two time periods:
  %    Pre SAL: T = 10 to 30 h (after RI, before encounter SAL)
  %        SAL: T = 40 to 60 h (during SAL)

  PreSalTstart = 10;
  PreSalTend   = 30;
  SalTstart    = 40;
  SalTend      = 60;

  % Description of cross sections
  XsectionList = {
    % in_file in_var pre_sal_out_var sal_out_var
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed_t_fa' '/ps_speed_t_fa' '/s_speed_t_fa' }
    { 'DIAGS/hist_meas_speed_<CASE>.h5' '/avg_speed_t_wm' '/ps_speed_t_wm' '/s_speed_t_wm' }
    };
  Nsets = length(XsectionList);


  for icase = 1:Ncases
    Case = Cases{icase};

    fprintf('**********************************************************\n');
    fprintf('Generating storm cross-sections for case: %s\n', Case);
    fprintf('\n');

    % Place all vars into one case specific output file.
    % If the file exists, remove it so that the HDF5 commands
    % can effectively re-create datasets.
    OutFile = sprintf('DIAGS/storm_xsections_%s.h5', Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    for iset = 1:Nsets
      InFile      = regexprep(XsectionList{iset}{1}, '<CASE>', Case);
      InVname     = XsectionList{iset}{2};
      PreSalVname = XsectionList{iset}{3};
      SalVname    = XsectionList{iset}{4};

      % Read in the variables, split into pre-SAL and SAL sections
      fprintf('  Reading: %s (%s)\n', InFile, InVname);

      VAR = squeeze(h5read(InFile, InVname));
      X   = squeeze(h5read(InFile, '/x_coords'));
      Y   = squeeze(h5read(InFile, '/y_coords'));
      Z   = squeeze(h5read(InFile, '/z_coords'));
      T   = squeeze(h5read(InFile, '/t_coords'));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % If this is the first set, write out the coordinates into the output file
      % so that subsequent variables can have these coordinates attached.
      if (iset == 1)
        Xname = '/x_coords';
        Yname = '/y_coords';
        Zname = '/z_coords';
        Tname = '/t_coords';
  
        CreateDimensionsXyzt(OutFile, X, Y, Z, [ 1 ], Xname, Yname, Zname, Tname);
        % Add COARDS annotations
        NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
      end

      % Determine the number of dimensions. When you have a vector (size is [ 1 n ] or [ n 1 ]),
      % Ndims will be set to 2 when you really want it to be 1.
      Ndims = ndims(VAR);
      if (Ndims == 2)
        if ((size(VAR,1) == 1) || (size(VAR,2) == 1))
          Ndims = 1
        end
      end
    
      % Form time indices for selecting the pre-SAL and SAL time periods
      SIM_TIME = (T ./ 3600) - 42;  % in hours, 0 through 60

      PT1 = find(SIM_TIME >= PreSalTstart, 1, 'first');
      PT2 = find(SIM_TIME <= PreSalTend,   1, 'last');

      ST1 = find(SIM_TIME >= SalTstart, 1, 'first');
      ST2 = find(SIM_TIME <= SalTend,   1, 'last');

      % Want to average across time, which is always the last dimension.
      % Reshape the output variables to 4D for Grads viewing.
      % Radius is the x dimension.
      switch(Ndims)
        case 1
          % var is 1D, (t)
          % results is a scalar
          P_VAR = squeeze(mean(VAR(PT1:PT2)));
          S_VAR = squeeze(mean(VAR(ST1:ST2)));

          DimOrder = { };
        case 2
          % var is 2D, (r,t)
          % result is (x)
          P_VAR = squeeze(mean(VAR(:,PT1:PT2),2));
          S_VAR = squeeze(mean(VAR(:,ST1:ST2),2));

          DimOrder = { 'x' };
        case 3
          % var is 3D, (r,z,t)
          % result is (x,z)
          P_VAR = squeeze(mean(VAR(:,:,PT1:PT2),3));
          S_VAR = squeeze(mean(VAR(:,:,ST1:ST2),3));

          DimOrder = { 'x' 'z' };
      end

      % Collect output vars into one case specific file.
      fprintf('  Writing: %s (%s, %s)\n', OutFile, PreSalVname, SalVname);
      h5create(OutFile, PreSalVname, size(P_VAR));
      h5create(OutFile, SalVname,    size(S_VAR));

      h5write(OutFile, PreSalVname, P_VAR);
      h5write(OutFile, SalVname,    S_VAR);

      AttachDimensionsXyzt(OutFile, PreSalVname, DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, SalVname,    DimOrder, Xname, Yname, Zname, Tname);

      fprintf('\n');
    end
  end
end
