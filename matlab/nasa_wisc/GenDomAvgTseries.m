function [ ] = GenDomAvgTseries()
% GenDomAvgTseries generate time series of domain averages

  Ddir = 'DIAGS';

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  % input_file input_dataset case_name outfile_prefix output_dataset
  VarSets = {
%    { 'HDF5/RCE70_OLD_RI/vint_vapor-a-AS-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE70_OLD_RI' 'avg_precip_water' 'precip_water' }
%    { 'HDF5/MATT/vint_vapor-a-AS-2012-01-01-000000-g1.h5'         'vertint_vapor' 'MATT'         'avg_precip_water' 'precip_water' }
%    { 'HDF5/RCE50_OLD_RI/vint_vapor-a-AS-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE50_OLD_RI' 'avg_precip_water' 'precip_water' }

%    { 'HDF5/RCE50_RECT/vint_vapor-a-AS-2012-01-01-000000-g1.h5'      'vertint_vapor' 'RCE50_RECT'      'avg_precip_water' 'precip_water' }
%    { 'HDF5/RCE50_RECT_S300/vint_vapor-a-AS-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE50_RECT_S300' 'avg_precip_water' 'precip_water' }
%    { 'HDF5/RCE50_RECT_S303/vint_vapor-a-AS-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE50_RECT_S303' 'avg_precip_water' 'precip_water' }
%    { 'HDF5/RCE50_SQ/vint_vapor-a-AS-2012-01-01-000000-g1.h5'        'vertint_vapor' 'RCE50_SQ'        'avg_precip_water' 'precip_water' }
%    { 'HDF5/RCE_MATT/vint_vapor-a-AS-2012-01-01-000000-g1.h5'        'vertint_vapor' 'RCE_MATT'        'avg_precip_water' 'precip_water' }
%    { 'HDF5/RCE_BASE/vint_vapor-a-AS-2012-01-01-000000-g1.h5'        'vertint_vapor' 'RCE_BASE'        'avg_precip_water' 'precip_water' }

%    { 'HDF5/RCE_CNTL/HDF5_NZ75/vint_vapor-RCE_CNTL-AC-2012-01-01-000000-g1.h5'   'vertint_vapor' 'RCE_CNTL_OLD_DIFF'        'avg_precip_water' 'precip_water' }
%    { 'HDF5/RCE_CNTL/HDF5/vint_vapor-RCE_CNTL-AC-2012-01-01-000000-g1.h5'        'vertint_vapor' 'RCE_CNTL_NEW_DIFF'        'avg_precip_water' 'precip_water' }

%    { 'HDF5/RCE_BASE/HDF5/vint_vapor-a-AC-2012-01-01-000000-g1.h5'                  'vertint_vapor' 'RCE_BASE'      'avg_precip_water' 'precip_water' }
%    { 'HDF5/RCE_EXP_S50LN/HDF5/vint_vapor-RCE_EXP_S50LN-AC-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE_EXP_S50LN' 'avg_precip_water' 'precip_water' }
%    { 'HDF5/RCE_EXP_S70LY/HDF5/vint_vapor-RCE_EXP_S70LY-AC-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE_EXP_S70LY' 'avg_precip_water' 'precip_water' }
%    { 'HDF5/RCE_EXP_S70LN/HDF5/vint_vapor-RCE_EXP_S70LN-AC-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE_EXP_S70LN' 'avg_precip_water' 'precip_water' }
%    { 'HDF5/RCE_EXP_S70MY/HDF5/vint_vapor-RCE_EXP_S70MY-AC-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE_EXP_S70MY' 'avg_precip_water' 'precip_water' }

%    { 'HDF5/RCE_EXP_S50LN_470/HDF5/vint_vapor-RCE_EXP_S50LN_470-AC-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE_EXP_S50LN_470' 'avg_precip_water' 'precip_water' }
%    { 'HDF5/RCE_EXP_S50LN_CG/HDF5/vint_vapor-RCE_EXP_S50LN_CG-AC-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE_EXP_S50LN_CG' 'avg_precip_water' 'precip_water' }
    { 'HDF5/RCE_EXP_S50LN_CGHZ/HDF5/vint_vapor-RCE_EXP_S50LN_CGHZ-AC-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE_EXP_S50LN_CGHZ' 'avg_precip_water' 'precip_water' }
    { 'HDF5/RCE_EXP_S50LN_THIN/HDF5/vint_vapor-RCE_EXP_S50LN_THIN-AC-2012-01-01-000000-g1.h5' 'vertint_vapor' 'RCE_EXP_S50LN_THIN' 'avg_precip_water' 'precip_water' }
    };
  Nset = length(VarSets);


  fprintf('***************************************************************\n');
  fprintf('Generating domain average time series:\n');

  for iset = 1:Nset
    InFile     = VarSets{iset}{1};
    InVname    = VarSets{iset}{2};
    Case       = VarSets{iset}{3};
    OutFprefix = VarSets{iset}{4};
    OutVname   = VarSets{iset}{5};

    OutFile = sprintf('%s/%s_%s.h5', Ddir, OutFprefix, Case);

    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Get entire field from input file
    % Do horizontal domain average
    %   if 2D field then result is a single number
    %   if 3D field then result is a vertical profile
    fprintf('    Reading: %s (%s)\n', InFile, InVname);
 
    HDATA = squeeze(hdf5read(InFile, InVname));
    Z     = squeeze(hdf5read(InFile, 'z_coords'));
    T     = squeeze(hdf5read(InFile, 't_coords'));
    Nz = length(Z);
    Nt = length(T);


    % Data is either (x,y,t) or (x,y,z,t) so calculate means on the first two dimensions and
    % write out the result
    AVG = squeeze(mean(HDATA,1));
    AVG = squeeze(mean(AVG,1));

    % Z will be the appropriate size
    OutVar = reshape(AVG, [ 1 1 Nz Nt ]);

    % Write out dummy coordinate values to keep ReadSelectXyzt happy
    Xdummy = 1;
    Ydummy = 1;

    fprintf('    Writing: %s (%s)\n', OutFile, OutVname);
    fprintf('\n');

    hdf5write(OutFile, OutVname, OutVar); 
    hdf5write(OutFile, '/x_coords', Xdummy, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', Ydummy, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', Z,      'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', T,      'WriteMode', 'append');
  end
end
