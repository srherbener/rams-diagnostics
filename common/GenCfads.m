function [ ] = GenCfads( ConfigFile )
% GenCfads generate CFADs from histogram data
%

% Read in the config data
[ Config ] = ReadConfig(ConfigFile);

CfadDir = Config.CfadDir;
ControlCase = Config.ControlCase;

% Write out the results into OutFile. Build the
% directory for OutFile in case it doesn't exist.
if (exist(CfadDir, 'dir') ~= 7)
    mkdir(CfadDir);
end

for icase = 1:length(Config.Cases)
  Case = Config.Cases(icase).Cname;
  for icfad = 1: length(Config.Cfads)
    Name = Config.Cfads(icfad).Name;
    InDir = Config.Cfads(icfad).InDir;
    Fprefix = Config.Cfads(icfad).Fprefix;
    Vname = Config.Cfads(icfad).Rvar;

    fprintf('***********************************************************************\n');
    fprintf('Generating CFAD:\n');
    fprintf('  Name: %s\n', Name);
    fprintf('  Variable: %s\n', Vname);
    fprintf('  Case: %s\n', Case);
    fprintf('\n');

    InFile  = sprintf('%s/%s_%s.h5', InDir, Fprefix, Case);
    OutFile = sprintf('%s/%s_%s.h5', CfadDir, Name, Case);
    Hdset   = sprintf('/%s', Vname);

    % Read in the histogram data. HIST will be organized as (r,b,z,t) where
    %    r --> radius
    %    b --> histogram bins
    %    z --> heights
    %    t --> time
    fprintf('Reading file: %s\n', InFile);
    fprintf('\n');
    R = hdf5read(InFile, '/x_coords');
    B = hdf5read(InFile, '/y_coords');
    Z = hdf5read(InFile, '/z_coords');
    T = hdf5read(InFile, '/t_coords');
    HIST = hdf5read(InFile, Hdset);
    [ Nr, Nb, Nz, Nt ] = size(HIST);

    % Change the bins from counts to frequency (fraction of sum of bins)
    %   1) Find the sums along the bins dimension, result is (r,1,z,t).
    %      In the case where all bin values are zero (the filtering mechanism
    %      caused no points to be selected), sum will be zero which will
    %      cause a divide by zero downstream. To prevent this, change all the
    %      zero sums to one. Since all of the bin values were zero in the first
    %      place the CFAD will end up with all zeros which is okay.
    %   2) Use repmat to tile the sums together to get back to the
    %      original size of the bins dimension where the elements
    %      for a given r,z,t are all the sum of that vector from the HIST array.
    %   3) Do an element wise divide to change each bin value to a frequency.
    SUMS = sum(HIST,2);
    SUMS(SUMS == 0) = 1;
    SUMS = repmat(SUMS, [ 1 Nb 1 1 ]);
    CFAD = HIST ./ SUMS;

    fprintf('Writing file: %s\n', OutFile);
    fprintf('\n');
    hdf5write(OutFile, Hdset, CFAD);
    hdf5write(OutFile, '/x_coords', R, 'WriteMode', 'append');
    hdf5write(OutFile, '/y_coords', B, 'WriteMode', 'append');
    hdf5write(OutFile, '/z_coords', Z, 'WriteMode', 'append');
    hdf5write(OutFile, '/t_coords', T, 'WriteMode', 'append');

    fprintf('\n');
  end
end

end
