function [ ] = GenMeanHydroParams()
% GenMeanHydroParams generate mean hydrometeor parameters for one-moment sim
%                    based on data from two-moment sim

  Ddir = 'DIAGS';

  % make sure output directory exists
  if (exist(Ddir, 'dir') ~= 7)
    mkdir(Ddir);
  end

  Case = 'RCE_S300';  % 2 moment case

  % mix_ratio_file mix_ratio_dataset number_conc_file number_conc_dataset
  VarSets = {
    { 1 'HDF5/TsAveragedData/hda_cloud_0p01_<CASE>.h5' '/hda_cloud' 'HDF5/TsAveragedData/hda_cloud_num_0p01_<CASE>.h5' '/hda_cloud_num' 'avg_param_cloud' '/cloud' }
    };
  Nset = length(VarSets);

  % reference density taken from RAMS output
  Dens = [
    1.152
    1.146
    1.140
    1.133
    1.126
    1.119
    1.111
    1.102
    1.093
    1.083
    1.072
    1.062
    1.050
    1.038
    1.026
    1.012
    0.998
    0.983
    0.969
    0.953
    0.937
    0.920
    0.903
    0.885
    0.867
    0.847
    0.828
    0.809
    0.788
    0.767
    0.745
    0.724
    0.700
    0.675
    0.651
    0.628
    0.603
    0.578
    0.553
    0.528
    0.503
    0.478
    0.453
    0.429
    0.406
    0.385
    0.364
    0.344
    0.326
    0.307
    0.290
    0.273
    0.257
    0.241
    0.226
    0.212
    0.198
    0.183
    0.164
    0.148
    0.133
    0.122
    0.112
    0.102
    0.094
    0.086
    0.079
    0.072
    0.066
    0.060
    0.055
    0.050
    0.046
    0.042
    0.039
    ];

  % RAMS settings for power law mass-diameter relationships for hydrometeors
  % 
  % ! shape      cfmas   pwmas      cfvt    pwvt     dmb0      dmb1
  % !---------------------------------------------------------------------- 
  %     .5,      524.,     3.,    3173.,     2.,   2.e-6,   50.e-6,  & !cloud
  %     .5,      524.,     3.,     149.,     .5,   .1e-3,    5.e-3,  & !rain
  %   .179,     110.8,   2.91,  5.769e5,   1.88,  15.e-6,  125.e-6,  & !pris col
  %   .179,  2.739e-3,   1.74,  188.146,   .933,   .1e-3,   10.e-3,  & !snow col
  %     .5,      .496,    2.4,    3.084,     .2,   .1e-3,   10.e-3,  & !aggreg
  %     .5,      157.,     3.,     93.3,     .5,   .1e-3,    5.e-3,  & !graup
  %     .5,      471.,     3.,     161.,     .5,   .8e-3,   10.e-3,  & !hail 
  %   .429,     .8854,    2.5,     316.,   1.01,      00,       00,  & !pris hex
  %  .3183,   .377e-2,     2.,     316.,   1.01,      00,       00,  & !pris den
  %  .1803,   1.23e-3,    1.8,  5.769e5,   1.88,      00,       00,  & !pris ndl
  %     .5,     .1001,  2.256,   3.19e4,   1.66,      00,       00,  & !pris ros
  %   .429,     .8854,    2.5,    4.836,    .25,      00,       00,  & !snow hex
  %  .3183,   .377e-2,     2.,    4.836,    .25,      00,       00,  & !snow den
  %  .1803,   1.23e-3,    1.8,  188.146,   .933,      00,       00,  & !snow ndl
  %     .5,     .1001,  2.256,  1348.38,  1.241,      00,       00,  & !snow ros
  %     .5,      524.,     3.,    3173.,     2.,  65.e-6,  100.e-6/    !drizzle
  %  
  % Use pris col and snow den habits
  %   1 cloud
  %   2 rain
  %   3 pris
  %   4 snow
  %   5 aggregates
  %   6 grapel
  %   7 hail
  %   8 drizzle
  Coeffs = [
    524.0
    524.0
    110.8
      0.377e-2
      0.496
    157.0
    471.0
      0.8854
    ];
  Powers = [
    3.0
    3.0
    2.91
    2.0
    2.4
    3.0
    3.0
    2.5
    ];

  fprintf('***************************************************************\n');
  fprintf('Generating parameter averages:\n');
  fprintf('  Case: %s\n', Case);
  fprintf('\n');


  for iset = 1:Nset
    Index      = VarSets{iset}{1};
    MassFile   = regexprep(VarSets{iset}{2}, '<CASE>', Case);
    MassVname  = VarSets{iset}{3};
    NumFile    = regexprep(VarSets{iset}{4}, '<CASE>', Case);
    NumVname   = VarSets{iset}{5};
    OutFile    = sprintf('%s/%s_%s.h5', Ddir, VarSets{iset}{6}, Case);
    OutVprefix = VarSets{iset}{7};

    OutDiamVname = sprintf('%s_diam', OutVprefix);
    OutNumVname  = sprintf('%s_num', OutVprefix);

    OutDiamProfVname = sprintf('%s_prof', OutDiamVname);
    OutNumProfVname = sprintf('%s_prof', OutNumVname);

    Cf  = Coeffs(Index);
    Pwr = Powers(Index);

    fprintf('  Hydrometeor Index: %d\n', Index);
    fprintf('    Coefficient: %f\n', Cf);
    fprintf('    Power: %f\n', Pwr);
    fprintf('\n');

    % Read in averages, format is (2,z,t) where
    %   (1,z,t) are the sums (per level)
    %   (2,z,t) are the counts (per level)
    fprintf('  Reading: %s (%s)\n', MassFile, MassVname);
    fprintf('  Reading: %s (%s)\n', NumFile, NumVname);
    fprintf('\n');

    MDATA = squeeze(h5read(MassFile, MassVname));
    NDATA = squeeze(h5read(NumFile, NumVname));
    Z     = squeeze(h5read(NumFile, '/z_coords'));

    MIXR  = squeeze(MDATA(1,:,:));  % sum of mass for each level and each time step
    CONC  = squeeze(NDATA(1,:,:));  % sum of number concentration for each level and each time step
    Nconc = squeeze(NDATA(2,:,:));  % count of selected grids that went into CONC

    MIXR = MIXR .* 1e-3; % convert from g/kg to kg/kg

    [ Nz, Nt ] = size(MIXR);

    % Get stretched vertical grid spacing
    DeltaZ = Z(2:end) - Z(1:end-1);
    DeltaZ = [ DeltaZ(1) DeltaZ' ]';  % just repeat the lower delta for the first entry

    % MIXR is kg/kg(air) and CONC is #/kg(air). The horizontal grid space is constant
    % everywhere, but the vertical grid space varies. Also the density varies with height.
    % In order to account for this create weights for averaging that are based on the
    % relative amounts of air at each vertical level. The volume changes by a factor of
    % DeltaZ(k)/DeltaZ(1) at each level, and the density changes by a factor of
    % Dens(k)/Dens(1) at each level. The overall factor is then:
    %
    %     Factor(k) = [ DeltaZ(k) * Dens(k) ] / [ DeltaZ(1) * Dens(1) ]
    %
    % Create WEIGHTS array with Factor as first dimension repeated along Time as
    % second dimension, that is (z,t). 
    Factor = (DeltaZ .* Dens) ./ (DeltaZ(1) * Dens(1));
    WEIGHTS = repmat(Factor, [ 1 Nt ]);

    % The volume of grid cells and the density of air are constant for a given
    % level, so don't apply the weights to get the profiles. But do apply weights
    % For total domain averages.

    % Assume power law for mass-diameter relationships. 
    %  m = Cf * D ^ Pwr
    %
    %  D = ( m / Cf ) ^ (1/Pwr)
    %    m = MIXR / CONC (average mass per particle, D will then be mean diameter)
    %
    % MIXR, CONC, and Nconc are (z,t)
    AVG_MASS = MIXR ./ CONC;

    AVG_DIAM_PROF = ( AVG_MASS ./ Cf) .^ (1 / Pwr);
    AVG_CONC_PROF = CONC ./ Nconc; 

    % Calculate average of domain mass by dividing the weighted sum of mixing ratio
    % over the domain by the weighted sum of number concentration over the domain
    % (as opposed to a weighted mean of mass calculated at every grid cell).
    AVG_DOM_MASS = squeeze(sum((MIXR .* WEIGHTS), 1)) ./ squeeze(sum((CONC .* WEIGHTS), 1));
 
    AVG_DIAM = (AVG_DOM_MASS ./ Cf) .^ (1 / Pwr);
    AVG_CONC  = squeeze(sum((CONC .* WEIGHTS), 1)) ./ squeeze(sum(WEIGHTS, 1));

    % Output
    fprintf('  Writing: %s (%s)\n', OutFile, OutDiamVname);
    fprintf('  Writing: %s (%s)\n', OutFile, OutNumVname);
    fprintf('  Writing: %s (%s)\n', OutFile, OutDiamProfVname);
    fprintf('  Writing: %s (%s)\n', OutFile, OutNumProfVname);
    fprintf('\n');

    % Remove existing file so that subsequent create/write command will not
    % have to deal with replacing data.
    if (exist(OutFile, 'file') == 2)
      delete(OutFile);
    end

    h5create(OutFile, OutDiamVname, size(AVG_DIAM));
    h5create(OutFile, OutNumVname, size(AVG_CONC));
    h5create(OutFile, OutDiamProfVname, size(AVG_DIAM_PROF));
    h5create(OutFile, OutNumProfVname, size(AVG_CONC_PROF));

    h5write(OutFile, OutDiamVname, AVG_DIAM);
    h5write(OutFile, OutNumVname, AVG_CONC);
    h5write(OutFile, OutDiamProfVname, AVG_DIAM_PROF);
    h5write(OutFile, OutNumProfVname, AVG_CONC_PROF);
  end
end
