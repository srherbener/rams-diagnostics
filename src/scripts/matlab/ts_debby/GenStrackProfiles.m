function [ ] = GenStrackProfiles()

  Ddir = 'DIAGS';
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  Cases = {
   'TSD_SAL_DUST'
   'TSD_NONSAL_DUST'
   };
  Ncases = length(Cases);

  ControlCase = 'TSD_SAL_DUST';

  % From TSD_SAL_DUST RAMS output (Reference density)
  RhoAir = [
    1.196
    1.191
    1.185
    1.179
    1.172
    1.165
    1.157
    1.147
    1.136
    1.124
    1.111
    1.097
    1.082
    1.066
    1.051
    1.034
    1.017
    1.000
    0.982
    0.963
    0.945
    0.925
    0.905
    0.883
    0.860
    0.833
    0.804
    0.773
    0.741
    0.711
    0.679
    0.647
    0.612
    0.577
    0.541
    0.505
    0.469
    0.432
    0.394
    0.355
    0.316
    0.279
    0.244
    0.210
    0.179
    0.150
    0.126
    0.105
    0.087
    0.073
    0.062
    0.052
    0.044
    0.038
    0.032
    0.027
    ];


  % Selection
  %
  % Along x-axis:
  CenterStart = 797; % 817 km is center of x range
  CenterEnd   = 837;

  % Along t-axis:
  PreSalStart = 10;
  PreSalEnd   = 30;

  SalStart = 40;
  SalEnd   = 60;

  % Description of profiles
  ProfileList = {
    % in_file in_var in_select out_var 
    %
    % data is (x,z,t)

    { 'XsectionData/strack_d1_num_<CASE>.h5' '/d1_num' 'center' '' '/center_d1_num' }
    { 'XsectionData/strack_d2_num_<CASE>.h5' '/d2_num' 'center' '' '/center_d2_num' }

    { 'XsectionData/strack_tracer1_<CASE>.h5' '/tracer1' 'center' '' '/center_tracer1' }
    { 'XsectionData/strack_tracer2_<CASE>.h5' '/tracer2' 'center' '' '/center_tracer2' }

    { 'XsectionData/strack_cloud_<CASE>.h5' '/cloud' 'center' '' '/center_cloud' }
    { 'XsectionData/strack_rain_<CASE>.h5'  '/rain'  'center' '' '/center_rain'  }

    };
  Nsets = length(ProfileList);


  for icase = 1:Ncases
    Case = Cases{icase};

    fprintf('**********************************************************\n');
    fprintf('Generating strack profiles for case: %s\n', Case);
    fprintf('\n');

    % Place all vars into one case specific output file.
    % If the file exists, remove it so that the HDF5 commands
    % can effectively re-create datasets.
    OutFile = sprintf('DIAGS/strack_profiles_%s.h5', Case);
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    for iset = 1:Nsets
      InFile   = regexprep(ProfileList{iset}{1}, '<CASE>', Case);
      InVname  = ProfileList{iset}{2};
      Xselect  = ProfileList{iset}{3};
      Tselect  = ProfileList{iset}{4};
      OutVname = ProfileList{iset}{5};

      ControlInFile = regexprep(ProfileList{iset}{1}, '<CASE>', ControlCase);

      OutDiffVname = sprintf('%s_diff', OutVname);

      % Read in the variables, VAR will be (x,z,t)
      fprintf('  Reading: %s (%s)\n', InFile, InVname);
      fprintf('  Reading: %s (%s)\n', ControlInFile, InVname);
      fprintf('\n');

      VAR = squeeze(h5read(InFile, InVname));
      X   = squeeze(h5read(InFile, '/x_coords'));
      Y   = squeeze(h5read(InFile, '/y_coords'));
      Z   = squeeze(h5read(InFile, '/z_coords'));
      T   = squeeze(h5read(InFile, '/t_coords'));

      CNTL_VAR = squeeze(h5read(ControlInFile, InVname));
      SIM_T = (T ./ 3600) - 42;

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % Apply selection and take average across x dimension
      % VAR will be (z,t) after this section

      switch(Xselect)
        case 'center'
          X1 = find(X >= CenterStart, 1, 'first');
          X2 = find(X <= CenterEnd,   1, 'last');
          ReduceX = true;

        otherwise
          X1 = 1;
          X2 = Nx;
          ReduceX = false;
      end

      switch(Tselect)
        case 'pre_sal'
          T1 = find(X >= PreSalStart, 1, 'first');
          T2 = find(X <= PreSalEnd,   1, 'last');
          ReduceT = true;

        case 'sal'
          T1 = find(X >= SalStart, 1, 'first');
          T2 = find(X <= SalEnd,   1, 'last');
          ReduceT = true;

        otherwise
          T1 = 1;
          T2 = Nt;
          ReduceT = false;
      end

      fprintf('  Reduction:\n');
      if (ReduceX)
        fprintf('    X range: %s -> From %.2f to %.2f (%d, %d)\n', Xselect, X(X1), X(X2), X1, X2);

        % (x,z,t) -> (z,t)
        VAR      = squeeze(nanmean(VAR(X1:X2,:,T1:T2), 1));
        CNTL_VAR = squeeze(nanmean(CNTL_VAR(X1:X2,:,T1:T2), 1));
        Vsize = [ Nz Nt ];
        DimOrder = { 'z' 't' };
      end
      if (ReduceT)
        fprintf('    T range: %s -> From %.2f to %.2f (%d, %d)\n', Tselect, SIM_T(T1), SIM_T(T2), T1, T2);
        if (ReduceX)
          % (z,t) -> (z)
          VAR      = squeeze(nanmean(VAR(:,T1:T2), 2));
          CNTL_VAR = squeeze(nanmean(CNTL_VAR(:,T1:T2), 2));
          Vsize = Nz;
          DimOrder = { 'z' };
        else
          % (x,z,t) -> (x,z)
          VAR      = squeeze(nanmean(VAR(X1:X2,:,T1:T2), 3));
          CNTL_VAR = squeeze(nanmean(CNTL_VAR(X1:X2,:,T1:T2), 3));
          Vsize = [ Nx Nz ];
          DimOrder = { 'x' 'z' };
        end
      end
      fprintf('\n');

      % Change nans to zeros to make the profiles look better on plots.
      VAR(isnan(VAR)) = 0;
      CNTL_VAR(isnan(CNTL_VAR)) = 0;

      DIFF_VAR = VAR - CNTL_VAR;

      % If this is the first set, write out the coordinates into the output file
      % so that subsequent variables can have these coordinates attached.
      if (iset == 1)
        Xname = '/x_coords';
        Yname = '/y_coords';
        Zname = '/z_coords';
        Tname = '/t_coords';
  
        CreateDimensionsXyzt(OutFile, X, Y, Z, T, Xname, Yname, Zname, Tname);
        % Add COARDS annotations
        NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);
      end

      % Collect output vars into one case specific file.
      fprintf('  Writing: %s (%s)\n', OutFile, OutVname);
      h5create(OutFile, OutVname, Vsize);
      h5write(OutFile, OutVname, VAR);

      fprintf('  Writing: %s (%s)\n', OutFile, OutDiffVname);
      h5create(OutFile, OutDiffVname, Vsize);
      h5write(OutFile, OutDiffVname, DIFF_VAR);

      AttachDimensionsXyzt(OutFile, OutVname, DimOrder, Xname, Yname, Zname, Tname);
      AttachDimensionsXyzt(OutFile, OutDiffVname, DimOrder, Xname, Yname, Zname, Tname);

      fprintf('\n');
    end
  end
end
