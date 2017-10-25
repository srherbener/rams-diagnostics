function [ ] = Hdf5ToBin (Fname, Vname, OutFile)
% Hdf5ToBin read in a variable from an hdf5 file and write it out in a binary format

  fprintf('Reading HDF5 file: %s, Dataset: %s\n', Fname, Vname);
  Var = hdf5read(Fname, Vname); 

  Dims = size(Var);
  Dnames = { 'x' 'y' 'z' 't' };
  fprintf('  Dimension sizes:\n');
  for i = 1:length(Dims)
    fprintf('    %s: %d\n', Dnames{i}, Dims(i));
  end
  fprintf('\n');

  fprintf('Writing binary data file: %s\n', OutFile);
  Fid = fopen(OutFile, 'w');
  fwrite(Fid, Var, 'real*4');
  fclose(Fid);

end
