function [ ] = ExtractHorizWinds()

  CaseList = {
    'TSD_SAL_DUST'
    'TSD_NONSAL_DUST'
    'TSD_SAL_NODUST'
    'TSD_NONSAL_NODUST'
    };
  Nc = length(CaseList);

  SampleList = {
    { 'HDF5/<CASE>/HDF5/u-<CASE>-AS-2006-08-20-120000-g3.h5' '/u' 10 3.5 '/u_t10_z35' }
    { 'HDF5/<CASE>/HDF5/u-<CASE>-AS-2006-08-20-120000-g3.h5' '/u' 20 3.5 '/u_t20_z35' }
    { 'HDF5/<CASE>/HDF5/u-<CASE>-AS-2006-08-20-120000-g3.h5' '/u' 30 3.5 '/u_t30_z35' }

    { 'HDF5/<CASE>/HDF5/v-<CASE>-AS-2006-08-20-120000-g3.h5' '/v' 10 3.5 '/v_t10_z35' }
    { 'HDF5/<CASE>/HDF5/v-<CASE>-AS-2006-08-20-120000-g3.h5' '/v' 20 3.5 '/v_t20_z35' }
    { 'HDF5/<CASE>/HDF5/v-<CASE>-AS-2006-08-20-120000-g3.h5' '/v' 30 3.5 '/v_t30_z35' }
    };
  Ns = length(SampleList);

  Xname = '/x_coords';
  Yname = '/y_coords';
  Zname = '/z_coords';
  Tname = '/t_coords';

  fprintf('Extracting horizontal wind samples:\n');
  for icase = 1:Nc
    Case = CaseList{icase};

    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Place all samples for one case in one file
    OutFname = sprintf('DIAGS/sample_hwind_%s.h5', Case);

    % Remove output file if it exists so that new dataset can be written.
    if (exist(OutFname, 'file') == 2)
      delete(OutFname);
    end

    for isamp = 1:Ns
      InFname    = regexprep(SampleList{isamp}{1}, '<CASE>', Case);
      InVname    = SampleList{isamp}{2};
      SampTime   = SampleList{isamp}{3};
      SampHeight = SampleList{isamp}{4};
      OutVname   = SampleList{isamp}{5};

      fprintf('    Reading: %s (%s)\n', InFname, InVname);
      X = squeeze(h5read(InFname, Xname));
      Y = squeeze(h5read(InFname, Yname));
      Z = squeeze(h5read(InFname, Zname));
      T = squeeze(h5read(InFname, Tname));

      Nx = length(X);
      Ny = length(Y);
      Nz = length(Z);
      Nt = length(T);

      % Find index corresponding to sample time, and use this to selectively read
      % just the sample out of the input file.
      SIM_T = T ./ 3600 - 42; % in hours, starting with zero
      Z_KM = Z ./ 1000;
      [ Dummy T1 ] = min(abs(SIM_T - SampTime));
      [ Dummy Z1 ] = min(abs(Z_KM - SampHeight));
      fprintf('    Sample time: %.1f h -> (%d, %.1f h)\n', SampTime, T1, SIM_T(T1));
      fprintf('    Sample height: %.1f km -> (%d, %.3f km)\n', SampHeight, Z1, Z_KM(Z1));

      Start = [ 1 1 Z1 T1 ];
      Count = [ Nx Ny 1 1 ];
      SAMP = squeeze(h5read(InFname, InVname, Start, Count));

      % Write out sample. If this is the first sample, then also write
      % out the coordinates.
      if (isamp == 1)
        CreateDimensionsXyzt(OutFname, X, Y, Z, T, Xname, Yname, Zname, Tname);
        NotateDimensionsXyzt(OutFname, Xname, Yname, Zname, Tname);
      end

      Vsize = [ Nx Ny ];
      DimOrder = { 'x' 'y' };

      fprintf('    Writing: %s (%s)\n', OutFname, OutVname);
      h5create(OutFname, OutVname, Vsize);
      h5write(OutFname, OutVname, SAMP);
      AttachDimensionsXyzt(OutFname, OutVname, DimOrder, Xname, Yname, Zname, Tname);
  
      fprintf('\n');
    end
  end
end
