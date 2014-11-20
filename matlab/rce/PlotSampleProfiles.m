function [ ] = PlotSampleProfiles(ConfigFile)
% PlotSampleProfiles plot vertical profile samples

  % Read the config file to get the structure of how the data is laid out in
  % the file system.
  [ Config ] = ReadConfig(ConfigFile);

  Ddir = Config.DiagDir;
  Pdir = Config.PlotDir;

  % make sure output directory exists
  if (exist(Pdir, 'dir') ~= 7)
    mkdir(Pdir);
  end

  % input_file input_dataset case_name outfile_prefix output_dataset
  PlotSets = {
      { 'RCE50_NEW_RI_MID' 'tempc' 'cloud' 'vapor' 240 'SZA 50, New Ri' }
      { 'RCE50_OLD_RI_MID' 'tempc' 'cloud' 'vapor' 240 'SZA 50, Old Ri' }
      { 'RCE70_OLD_RI_MID' 'tempc' 'cloud' 'vapor' 240 'SZA 70, Old Ri' }
    };
  Nset = length(PlotSets);


  fprintf('***************************************************************\n');
  fprintf('Generating vertical profile plots:\n');

  for iset = 1:Nset
    Case   = PlotSets{iset}{1};

    TempVar  = PlotSets{iset}{2};
    CloudVar = PlotSets{iset}{3};
    VaporVar = PlotSets{iset}{4};

    SelTime = PlotSets{iset}{5};

    Ptitle = PlotSets{iset}{6};

    InFile = sprintf('%s/SampleProfiles_%s.h5', Ddir, Case);

    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    % Read in TEMP, CLOUD, and VAPOR data
    % After squeezing, vars will be: (z,t)
    fprintf('    Reading: %s (%s, %s, %s)\n', InFile, TempVar, CloudVar, VaporVar);
    fprintf('      Selected time: %f.2\n', SelTime);
    fprintf('\n');

    TEMP = squeeze(hdf5read(InFile, TempVar));
    CLOUD = squeeze(hdf5read(InFile, CloudVar));
    VAPOR = squeeze(hdf5read(InFile, VaporVar));

    Z = squeeze(hdf5read(InFile, 'z_coords')) ./ 1000;  % km
    T = squeeze(hdf5read(InFile, 't_coords')) ./ 3600;  % h

    [ Nz Nt ] = size(TEMP);

    % Trim out the selected dimensions.
    Z1 = 2;
    Z2 = Nz;

    T1 = find(T >= SelTime, 1, 'first');

    TEMP = squeeze(TEMP(Z1:Z2,T1));
    CLOUD = squeeze(CLOUD(Z1:Z2,T1));
    VAPOR = squeeze(VAPOR(Z1:Z2,T1));

    Z = Z(Z1:Z2);
    
    % Create the plot
    PDATA = [ TEMP CLOUD VAPOR ];
    LegendText = { 'T (deg C)' 'Qc (g/kg)' 'Qv (g/kg)' };

    Ptitle = sprintf('%s, t = %d h', Ptitle, SelTime);

    Fig = figure;

    plot(PDATA,Z, 'LineWidth', 1.5);
    legend(LegendText);
    set(gca, 'FontSize', 20);
    title(Ptitle);
    ylabel('Height (km)');

    OutFile = sprintf('%s/SampleSounding_%s_T%d.png', Pdir, Case, SelTime);
    fprintf('    Writing: %s\n', OutFile);
    fprintf('\n');

    saveas(Fig, OutFile);
    close(Fig);
  end
end
