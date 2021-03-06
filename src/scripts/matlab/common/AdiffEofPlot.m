function [ ] = AdiffEofPlot(InFile, EofOutFile, PcOutFile, EsOutFile, Vname, Vunits, Ptitle, EofNum, NumEv, SelectData, Clevs, Cbounds)
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
Z = hdf5read(InFile, 'Height') / 1000; % convert to km
T = hdf5read(InFile, 'Time') / 3600; % convert to hr

% Translate the coordinate values in SelectData to their corresponding
% indices.
SelectIndices(1) = find(R >= SelectData(1), 1);
SelectIndices(2) = find(R <= SelectData(2), 1, 'last');
SelectIndices(3) = find(Z >= SelectData(3), 1);
SelectIndices(4) = find(Z <= SelectData(4), 1, 'last');

% Pick out the EOF/PC pair of interest
EOF = ALL_EOF(:,EofNum);
PC = ALL_PC(:,EofNum);

% Generate the plot
E_Title = Ptitle;
E_Xlabel = 'Radius (km)';
E_Ylabel = 'Height (km)';

P_Title = Ptitle;
P_Xlabel = 'Simulation Time (h)';
P_Ylabel = sprintf('%s (%s)', Vname, Vunits);

fprintf('Writing EOF plot to file: %s\n', EofOutFile);
Plot2dEofPc(EOF, PC, R, Z, Clevs, Cbounds, E_Title, E_Xlabel, E_Ylabel, P_Title, P_Xlabel, P_Ylabel, SelectIndices, EofOutFile, PcOutFile);

EsTitle = sprintf('First %d Eigenvalues of %s Difference: %s', NumEv, Vname, regexprep(Ptitle, '^PANEL:', ''));
fprintf('Writing eigenvalue spectrum plot to file: %s\n', EsOutFile);
PlotEigSpectrum( VarExpl, Err, NumEv, EsTitle, EsOutFile );

end
