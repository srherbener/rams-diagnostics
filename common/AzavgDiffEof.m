function [ ] = AzavgDiffEof( InFile1, InFile2, Vname, OutFile, DataSelect, UndefVal, Nstar )
% AzavgDiffEof generate the EOF/PCs and the Eigenvalue spectrum for a given dataset
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

% Read in the HDF5 data. The data will be organized as (x,y,z,t) after
% being read in.
fprintf('Reading %s from HDF5 file: %s\n', Vname, InFile1);
H5_Var1  = hdf5read(InFile1, Vname);
fprintf('Reading %s from HDF5 file: %s\n', Vname, InFile2);
H5_Var2  = hdf5read(InFile2, Vname);
fprintf('\n');

% convert the undefined value (-999) to nan
H5_Var1(H5_Var1 == UndefVal) = nan;
H5_Var2(H5_Var2 == UndefVal) = nan;

% Read in coordinate values from InFile1 (doesn't matter which one). Don't
% need y_coords since this is a dummy dimension.
F1_Xcoords = hdf5read(InFile1, 'x_coords');
F1_Zcoords = hdf5read(InFile1, 'z_coords');
F1_Tcoords = hdf5read(InFile1, 't_coords');

F2_Xcoords = hdf5read(InFile2, 'x_coords');
F2_Zcoords = hdf5read(InFile2, 'z_coords');
F2_Tcoords = hdf5read(InFile2, 't_coords');

% Grab the selection criteria out of DataSelect
F1X1 = find(F1_Xcoords >= DataSelect(1), 1);
F1X2 = find(F1_Xcoords <= DataSelect(2), 1, 'last');
F1Y1 = 1;
F1Y2 = 1;
F1Z1 = find(F1_Zcoords >= DataSelect(3), 1);
F1Z2 = find(F1_Zcoords <= DataSelect(4), 1, 'last');
F1T1 = find(F1_Tcoords >= DataSelect(5), 1);
F1T2 = find(F1_Tcoords <= DataSelect(6), 1, 'last');

F2X1 = find(F2_Xcoords >= DataSelect(1), 1);
F2X2 = find(F2_Xcoords <= DataSelect(2), 1, 'last');
F2Y1 = 1;
F2Y2 = 1;
F2Z1 = find(F2_Zcoords >= DataSelect(3), 1);
F2Z2 = find(F2_Zcoords <= DataSelect(4), 1, 'last');
F2T1 = find(F2_Tcoords >= DataSelect(5), 1);
F2T2 = find(F2_Tcoords <= DataSelect(6), 1, 'last');

% Convert to the 2D format that EOF analysis wants --> (time, obs)
Var1 = Xyzt2EofArray(H5_Var1, F1X1, F1X2, F1Y1, F1Y2, F1Z1, F1Z2, F1T1, F1T2);
Var2 = Xyzt2EofArray(H5_Var2, F2X1, F2X2, F2Y1, F2Y2, F2Z1, F2Z2, F2T1, F2T2);

% Create the difference values
% Wait until here to create the differences since the control
% run will typically have more time steps in it (making H5_Var1
% and H5_Var2 different sizes), and the data selection in
% Xyzt2EofArray will trim out equal sizes.
Var = Var1 - Var2;
%Var(isnan(Var)) = 0;

% Coordinate values
% Use the coordinates from the first file (not the control)
Radius = F1_Xcoords(F1X1:F1X2);
Height = F1_Zcoords(F1Z1:F1Z2);
Time   = F1_Tcoords(F1T1:F1T2);

% Run EOF and generate the eigenvalue spectrum
fprintf('Running EOF analysis: %s %s\n', Vname, InFile1);
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

