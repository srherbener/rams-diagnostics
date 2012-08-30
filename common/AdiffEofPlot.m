function [ ] = AdiffEofPlot(InFile, EofOutFile, EsOutFile, Vname, Vunits, Pcase, EofNum, NumEv, SelectData, Clevs, Cbounds)
% AdiffEofPlot plot an EOF/PC pair that was creaged by Azavg difference EOF routine.
%
% This function will take the EOF/PC data from Infile and plot the EOF/PC pair
% corresponding to EofNum.

% Read in the HDF5 data
fprintf('Reading EOF data from HDF5 file: %s\n', InFile);
ALL_EOF = hdf5read(InFile, 'EOF');
ALL_PC  = hdf5read(InFile, 'PC');
VarExpl = hdf5read(InFile, 'VarExpl');
Err     = hdf5read(InFile, 'Err');

R = hdf5read(InFile, 'Radius') / 1000; % convert to km
Z = hdf5read(InFile, 'Height');
T = hdf5read(InFile, 'Time') / 3600; % convert to hr

% Pick out the EOF/PC pair of interest
EOF = ALL_EOF(:,EofNum);
PC = ALL_PC(:,EofNum);

% Generate the plot
E_Title = sprintf('%s Difference EOF%d: %s', Vname, EofNum, Pcase);
E_Xlabel = 'Radius (km)';
E_Ylabel = 'Height (m)';

P_Title = sprintf('%s Difference PC%d: %s', Vname, EofNum, Pcase);
P_Xlabel = 'Time (hr)';
P_Ylabel = sprintf('%s (%s)', Vname, Vunits);

fprintf('Writing EOF plot to file: %s\n', EofOutFile);
Plot2dEofPc(EOF, PC, R, Z, Clevs, Cbounds, E_Title, E_Xlabel, E_Ylabel, P_Title, P_Xlabel, P_Ylabel, EofOutFile);

EsTitle = sprintf('First %d Eigenvalues of %s Difference: %s', NumEv, Vname, Pcase);
fprintf('Writing eigenvalue spectrum plot to file: %s\n', EsOutFile);
PlotEigSpectrum( VarExpl, Err, NumEv, EsTitle, EsOutFile );

end
