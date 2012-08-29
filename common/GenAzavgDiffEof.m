function [ ] = GenAzavgDiffEof( InFile1, Vname1, InFile2, Vname2, OutFile, DataSelect, UndefVal, Nstar )
% GenAzavgDiffEof generate the EOF/PCs and the Eigenvalue spectrum for a given dataset
%
%   This function will read in data contained in InFile1:Vname1 and InFile2:Vname2, select
%   a subset of the data according to the specs in 'DataSelect', create a difference
%   of those subsets (Var1 - Var2) and perform an EOF analysis on the difference
%   data set. The results will be placed in 'OutFile'. All files are HDF5 format
%   that comes out of revu.
%
%   It is assumed that the input data was created using the azavg program. This means
%   that the input data is organized as (x,y,z,t) with y being a dummy dimension (size
%   eqauls one) to keep GRADS happy.
%
%   DataSelect is a vector, (x1, x2, z1, z3, t1, t2), showing the ranges of x, z, t
%   to use for selection. Note that only a rectangular volume can be selected using
%   this scheme - hopefully, that's sufficient.
%
%   UndefVal is the numeric value representing an undefined data value.
%
%   Nstar is the estimated sample size for the eigenvalue spectrum calculation.
%

% Read in the HDF5 data. The data will be organized as (x,y,z,t) after being read in.
fprintf('Reading %s from HDF5 file: %s\n', Vname1, InFile1);
H5_Var1  = hdf5read(InFile1, Vname1);
fprintf('Reading %s from HDF5 file: %s\n', Vname2, InFile2);
H5_Var2  = hdf5read(InFile2, Vname2);
fprintf('\n');

% Create the difference values
H5_Var = H5_Var1 - H5_Var2;

% Read in coordinate values from InFile1 (doesn't matter which one). Don't
% need y_coords since this is a dummy dimension.
H5_Xcoords = hdf5read(InFile1, 'x_coords');
H5_Zcoords = hdf5read(InFile1, 'z_coords');

% Grab the selection criteria out of DataSelect
X1 = DataSelect(1);
X2 = DataSelect(2);
Z1 = DataSelect(3);
Z2 = DataSelect(4);
T1 = DataSelect(5);
T2 = DataSelect(6);

% Convert to the 2D format that EOF analysis wants --> (time, obs)
Var = Xyzt2EofArray(H5_Var, X1, X2, 1, 1, Z1, Z2, T1, T2);

% Coordinate values
Radius = H5_Xcoords(X1:X2);
Height = H5_Zcoords(Z1:Z2);

% convert the undefined value (-999) to nan
Var(Var == UndefVal) = nan;

% Run EOF and generate the eigenvalue spectrum
fprintf('Running EOF analysis:\n', Vname1, InFile1);
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

end

