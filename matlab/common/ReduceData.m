function [ ] = ReduceData(InFile,OutFile,Var,Xinc,Yinc,Zinc,Tinc)
% ReduceData - Read in an HDF5 file, Reduce Var, and write out the results
% 
% This routine assumes that Var is a 4D variable (x,y,z,t) and that the
% coordinates are in datasets called "/[xyzt]_coords".
%
% Reduction is just taking every nth point according to [XYZT]inc.
%

  fprintf('*****************************************************************\n');
  fprintf('Reducing HDF5 data:\n');
  fprintf('  Reading: %s (%s)\n', InFile, Var);
  fprintf('  Writing: %s (%s)\n', OutFile, Var);
  fprintf('\n');


  % If OutFile exists, clobber it
  if (exist(OutFile, 'file') == 2)
    delete(OutFile);
  end
  
  % First read in the coordinate data, reduce these according to the
  % increment specs.
  X = squeeze(h5read(InFile, '/x_coords'));
  Y = squeeze(h5read(InFile, '/y_coords'));
  Z = squeeze(h5read(InFile, '/z_coords'));
  T = squeeze(h5read(InFile, '/t_coords'));

  Nx = length(X);
  Ny = length(Y);
  Nz = length(Z);
  Nt = length(T);

  X_OUT = X(1:Xinc:end);
  Y_OUT = Y(1:Yinc:end);
  Z_OUT = Z(1:Zinc:end);
  T_OUT = T(1:Tinc:end);

  NxOut = length(X_OUT);
  NyOut = length(Y_OUT);
  NzOut = length(Z_OUT);
  NtOut = length(T_OUT);

  Xname = '/x_coords';
  Yname = '/y_coords';
  Zname = '/z_coords';
  Tname = '/t_coords';

  CreateDimensionsXyzt(OutFile, X_OUT, Y_OUT, Z_OUT, T_OUT, Xname, Yname, Zname, Tname);
  % Add COARDS annotations
  NotateDimensionsXyzt(OutFile, Xname, Yname, Zname, Tname);

  % Expect Var to be a large dataset. Read one time step at a time to keep
  % memory usage reasonable.

  h5create(OutFile, Var, [ NxOut, NyOut, NzOut, NtOut ]);

  count = 0;
  for it = 1:Tinc:Nt
    % VAR will be 3D (x,y,z) after this read
    VAR = squeeze(h5read(InFile, Var, [ 1 1 1 it ], [ Nx Ny Nz 1 ]));

    VAR_OUT = VAR(1:Xinc:end,1:Yinc:end,1:Zinc:end);

    count = count + 1;
    h5write(OutFile, Var, VAR_OUT, [ 1 1 1 count ], [ NxOut, NyOut, NzOut, 1 ]);

    if (mod(it, 10) == 0)
      fprintf('  Working... time step: %d\n', it);
    end   
  end
  fprintf('\n');

  fprintf('  Processed %d time steps\n', it-1);
  fprintf('\n');

  % Attach dimensions
  DimOrder = { 'x' 'y' 'z' 't' };
  AttachDimensionsXyzt(OutFile, Var, DimOrder, Xname, Yname, Zname, Tname);

end % function
