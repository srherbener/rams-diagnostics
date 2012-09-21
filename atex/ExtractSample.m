function [ ] = ExtractSample(ConfigFile)
% Script to extract sample data for a 4 panel plot
%

Config = ReadConfig(ConfigFile);

Tdir = Config.TsavgDir;
Ddir = Config.DiagDir;

SampleDir = 'SampleHdf5';

FileSuffix = '-AS-1999-02-10-040000-g1.h5';


% Read in the PCPRR data
%
% After reading in, the dimensions will be (x,y,t)

CCN = [ 50 50 1600 1600 ]'; % keep these the same length
SST = [ 298 303 298 303 ]';

for i = 1:length(CCN)
 Exps{i} = sprintf('z.atex.ccn%04d.sst%03d.gcn10m5', CCN(i), SST(i));
end

Vars = { 'pcprr';
'vint_cond';
'ctop_tempc' };

Rvars = { 'pcprr';
'vertint_cond';
'cloudtop_TS' };

% For this sim: 1999, Feb 10, 2300Z is time step 229
Tstep = 229;
Tstring = '1999-02-10T23:00:00Z';

% Write the time step into the file
h5_fout = sprintf('%s/SampleData.h5', Ddir);
hdf5write(h5_fout, '/Tstep', Tstep);
hdf5write(h5_fout, '/Time',Tstring, 'WriteMode', 'append');

% Copy the lat and lon values into the output file
h5_grid_fin = sprintf('%s/%s/%s-%s%s', SampleDir, Exps{1}, Vars{1}, Exps{1}, FileSuffix);
Lon = hdf5read(h5_grid_fin, '/x_coords');
Lat = hdf5read(h5_grid_fin, '/y_coords');
hdf5write(h5_fout, '/Lon', Lon, 'WriteMode', 'append');
hdf5write(h5_fout, '/Lat', Lat, 'WriteMode', 'append');

% Write out the CCN and SST values
hdf5write(h5_fout, '/CcnConcen', CCN, 'WriteMode', 'append');
hdf5write(h5_fout, '/Sst', SST, 'WriteMode', 'append');

% Rain rate is in the REVU var PCPRR
for i = 1:length(Exps)
  Exp = Exps{i};
  for j = 1:length(Vars)
    Var = Vars{j};
    h5_fin = sprintf('%s/%s/%s-%s%s', SampleDir, Exps{i}, Vars{j}, Exps{i}, FileSuffix);
    h5_dset_in = sprintf('/%s',Rvars{j});
    h5_dset_out = sprintf('/CCN_%d/SST_%d/%s',CCN(i),SST(i),Var);
    fprintf('Extracting Data:\n');
    fprintf('  Experiment: %s\n', Exp);
    fprintf('  Variable: %s\n', Var);
    fprintf('  Input file: %s\n',h5_fin);
    fprintf('\n');

    % Read the variable out of the file and pull out
    % the field corresponding to Tstep and write into
    % the output file.
    h5_data = hdf5read(h5_fin, h5_dset_in);

    Vdata = h5_data(:,:,Tstep);
    hdf5write(h5_fout, h5_dset_out, Vdata, 'WriteMode', 'append');
  end
end
