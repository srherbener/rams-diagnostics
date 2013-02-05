function [ ] = WriteData (Fname, Var)
% WriteData write array data out into a binary file

  Fid = fopen(Fname, 'w');
  fwrite(Fid, Var, 'real*4');
  fclose(Fid);

end
