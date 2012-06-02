% Script to extract sample data for a 4 panel plot
%

clear;

% Read in the PCPRR data
%
% After reading in, the dimensions will be (x,y,t)

Exps = { 'z.atex250m.100km.ccn0050.sst298';
'z.atex250m.100km.ccn0050.sst303';
'z.atex250m.100km.ccn1600.sst298';
'z.atex250m.100km.ccn1600.sst303' };

CCN = [ 50 50 1600 1600 ]'; % keep these in sync with order in Exps
SST = [ 298 303 298 303 ]';

Vars = { 'PCPRR';
'VERTINT_COND';
'CLOUDTOP_TEMPC' };

% For this sim: 1999, Feb 10, 2300Z is time step 229
Tstep = 229;
Tstring = '1999-02-10T23:00:00Z';

% Write the time step into the file
h5_fout = sprintf('DIAG/SampleData.h5');
hdf5write(h5_fout, '/Tstep', Tstep);
hdf5write(h5_fout, '/Time',Tstring, 'WriteMode', 'append');

% Rain rate is in the REVU var PCPRR
for i = 1:size(Exps,1)
  Exp = char(Exps(i));
  for j = 1:size(Vars,1)
    Var = char(Vars(j));
    h5_fin = sprintf('REVU/%s/%s.h5',Exp,Var);
    h5_dset_in = sprintf('/%s',Var);
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
