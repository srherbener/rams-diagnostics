function [ ] = GenAzavgDiffEof( InFile1, InFile2, Vname, OutFile, DataSelect, UndefVal, Nstar )
% GenAzavgDiffEof generate the EOF/PCs and the Eigenvalue spectrum for a given dataset
%
%   This function will read in data contained in InFile1:Vname and InFile2:Vname, select
%   a subset of the data according to the specs in 'DataSelect', create a difference
%   of those subsets (Var1 - Var2) and perform an EOF analysis on the difference
%   data set. The results will be placed in 'OutFile'. All files are HDF5 format
%   that comes out of revu.
%
%   It is assumed that the input data was created using the azavg program. This means
%   that the input data is organized as (x,y,z,t) with y being a dummy dimension (size
%   eqauls one) to keep GRADS happy.
%
%   DataSelect is a vector, (x1, x2, z1, z2, t1, t2), showing the ranges of x, z, t
%   to use for selection. Note that only a rectangular volume can be selected using
%   this scheme - hopefully, that's sufficient.
%
%   UndefVal is the numeric value representing an undefined data value.
%
%   Nstar is the estimated sample size for the eigenvalue spectrum calculation.
%

% Read in the HDF5 data. The data will be organized as (x,y,z,t) after being read in.
if ((strcmp(Vname, 'w_up')) || (strcmp(Vname, 'w_dn')))
  H5_Var = 'w';
else
  H5_Var = Vname;
end
fprintf('Reading %s from HDF5 file: %s\n', H5_Var, InFile1);
H5_Var1  = hdf5read(InFile1, H5_Var);
fprintf('Reading %s from HDF5 file: %s\n', H5_Var, InFile2);
H5_Var2  = hdf5read(InFile2, H5_Var);
fprintf('\n');

% convert the undefined value (-999) to nan
H5_Var1(H5_Var1 == UndefVal) = nan;
H5_Var2(H5_Var2 == UndefVal) = nan;

% Create the difference values
H5_Var = H5_Var1 - H5_Var2;

% Read in coordinate values from InFile1 (doesn't matter which one). Don't
% need y_coords since this is a dummy dimension.
H5_Xcoords = hdf5read(InFile1, 'x_coords');
H5_Zcoords = hdf5read(InFile1, 'z_coords');
H5_Tcoords = hdf5read(InFile1, 't_coords');

% Grab the selection criteria out of DataSelect
X1 = find(H5_Xcoords >= DataSelect(1), 1);
X2 = find(H5_Xcoords <= DataSelect(2), 1, 'last');
Y1 = 1;
Y2 = 1;
Z1 = find(H5_Zcoords >= DataSelect(3), 1);
Z2 = find(H5_Zcoords <= DataSelect(4), 1, 'last');
T1 = find(H5_Tcoords >= DataSelect(5), 1);
T2 = find(H5_Tcoords <= DataSelect(6), 1, 'last');

% Convert to the 2D format that EOF analysis wants --> (time, obs)
Var = Xyzt2EofArray(H5_Var, X1, X2, Y1, Y2, Z1, Z2, T1, T2);

% Coordinate values
Radius = H5_Xcoords(X1:X2);
Height = H5_Zcoords(Z1:Z2);
Time   = H5_Tcoords(T1:T2);

% Run EOF and generate the eigenvalue spectrum
fprintf('Running EOF analysis:\n', Vname, InFile1);
[ EOF, PC, EigVals ] = EofNan(Var,1);
[ Lambda, VarExpl, Err ] = EigenSpectrum(EigVals, Nstar);
fprintf('\n');

% Write out the results into OutFile. Build the
% directory for OutFile in case it doesn't exist.
[ Fdir, Fname, Fext ] = fileparts(OutFile);
if (exist(Fdir, 'dir') ~= 7)
    mkdir(Fdir);
end

fprintf('Saving EOF results in: %s\n', OutFile);
hdf5write(OutFile, '/EOF',     EOF);
hdf5write(OutFile, '/PC',      PC,      'WriteMode', 'append');
hdf5write(OutFile, '/Lambda',  Lambda,  'WriteMode', 'append');
hdf5write(OutFile, '/VarExpl', VarExpl, 'WriteMode', 'append');
hdf5write(OutFile, '/Err',     Err,     'WriteMode', 'append');
hdf5write(OutFile, '/Radius',  Radius,  'WriteMode', 'append');
hdf5write(OutFile, '/Height',  Height,  'WriteMode', 'append');
hdf5write(OutFile, '/Time',    Time,    'WriteMode', 'append');

end

