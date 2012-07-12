% Script to generate time domain averaged profiles.
%
% This script assumes that the time steps go in order when you
% concatenate together the results from each time directory
% in the order specified in Tdirs.
%
% This allows this script to build an array with each profile
% in a row and then call the mean() function to generated the
% overall avg profile.

clear;

Var = 'lh_vapt';
Hvar = 'latheatvapt';
OutDir = 'DIAGS';

% Create the output directory if it doesn't exist
if (exist(OutDir,'dir') == 0)
  mkdir(OutDir);
end

[ Exps, Tdirs ] = ReadConfig('Config');

fprintf('Running time, domain averaging on: %s\n', Var);
for i_exp = 1:length(Exps)
  fprintf('  Experiment: %s\n', char(Exps(i_exp)));

  for i_tdir = 1:length(Tdirs)
    % Form the directory: <Exp>/HDF5/<Tdir>
    Hdir = sprintf('%s/HDF5/%s', char(Exps(i_exp)), char(Tdirs(i_tdir)));
    Hfpat = sprintf('%s/%s-*', Hdir, Var);
    Hfile = dir(Hfpat);
    Hfile = Hfile.name;
    
    % Need to pre-pend the directory since dir returns the basename of the file
    Hfile = sprintf('%s/%s', Hdir, Hfile);
    fprintf('    Reading HDF5 file: %s\n', Hfile);

    % Read in data for the variable from this time directory
    Vdata = hdf5read(Hfile, Hvar);

    % Find the average for this piece
    [ AvgProf ] = DomainTimeAvgProfile(Vdata);

    % Put this into a 2D array, each row is a profile from each time directory
    if (i_tdir == 1)
      AllProfs = zeros(length(Tdirs), length(AvgProf)); % initialize the 2D array
    end
    AllProfs(i_tdir,:) = [ AvgProf ];
  end
  fprintf('\n');

  % Average all of the time directory
  TotAvgProf = mean(AllProfs,1);

  % Save the profile in an hdf5 file
  Hfile = sprintf('./%s/%s_tdavg_%s.h5', OutDir, Var, char(Exps(i_exp)));
  Hdset = sprintf('/%s',Var);
  fprintf('    Writing HDF5 file: %s\n', Hfile);
  hdf5write(Hfile, Hdset, TotAvgProf);
  fprintf('\n');
end
