function [ Cdata ] = ReadConfig ( Cfile )
% ReadConfig read diagnostic configuration file for a set of simulations
%
% This function will configuration file for processing
% a set of experimental simulation runs.
%
% The config file has the format:
%    Case: <directory>
%    TimeDir: <directory>
%
% Case: and TimeDir: define all of the directories where hdf5 files live
% 
% Case: defines separate experiments (simulations)
%
% TimeDir: defines multiple directories under a single experiment where
%       hdf5 files live
%
% If you have
%   Case: RUN1
%   Case: RUN2
%   Case: RUN3
%   TimeDir: T1
%   TimeDir: T2
%
% then the directories where the hdf5 files live are:
%
%   For RUN1:
%     RUN1/HDF5/T1
%     RUN1/HDF5/T2
%
%   For RUN2:
%     RUN2/HDF5/T1
%     RUN2/HDF5/T2
%
%   For RUN3:
%     RUN3/HDF5/T1
%     RUN3/HDF5/T2
%

% Read in entire file and place into a character array.
fid = fopen(char(Cfile));
InFile = textscan(fid, '%s', 'Delimiter', '\n'); % read in the entire line
fclose(fid);

% Convert to char array, entire file is in first element of InFile.
InLines = char(InFile{1});

% Grab keywords and their data and load into structs
i_case = 0;
i_tdir = 0;
i_aeof = 0;
i_aeofplot = 0;
i_tsplot = 0;
for i = 1:length(InLines)
  % Convert a line to a list of space separated fields
  [ Fields ] = Line2Fields(InLines(i,:), ' ');


  switch Fields{1}
    case 'AzavgDir:'
      Cdata.AzavgDir = Fields{2};
    case 'TsavgDir:'
      Cdata.TsavgDir = Fields{2};
    case 'PlotDir:'
      Cdata.PlotDir = Fields{2};
    case 'EofDir:'
      Cdata.EofDir = Fields{2};
    case 'ControlCase:'
      Cdata.ControlCase = Fields{2};
    case 'UndefVal:'
      Cdata.UndefVal = sscanf(Fields{2}, '%f');
    case 'Case:'
      i_case = i_case + 1;
      Cdata.Cases(i_case).Cname = Fields{2};
      Cdata.Cases(i_case).Pname = Fields{3};
    case 'TimeDir:'
      i_tdir = i_tdir + 1;
      Cdata.Tdirs{i_tdir} = Fields{2};
    case 'AzavgEofConfig:'
      Cdata.AzavgEofConfig.NumEv    = sscanf(Fields{2}, '%d');
      Cdata.AzavgEofConfig.Nstar    = sscanf(Fields{3}, '%d');
    case 'AzavgEof:'
      i_aeof = i_aeof + 1;
      Cdata.AzavgEof(i_aeof).Var  = Fields{2};
      Cdata.AzavgEof(i_aeof).Rmin = sscanf(Fields{3}, '%f');
      Cdata.AzavgEof(i_aeof).Rmax = sscanf(Fields{4}, '%f');
      Cdata.AzavgEof(i_aeof).Zmin = sscanf(Fields{5}, '%f');
      Cdata.AzavgEof(i_aeof).Zmax = sscanf(Fields{6}, '%f');
      Cdata.AzavgEof(i_aeof).Tmin = sscanf(Fields{7}, '%f');
      Cdata.AzavgEof(i_aeof).Tmax = sscanf(Fields{8}, '%f');
    case 'PlotExp:'
      Cdata.Pexp.Ename  = Fields{2};
      Cdata.Pexp.Tstart = sscanf(Fields{3}, '%d');
      Cdata.Pexp.Tend   = sscanf(Fields{4}, '%d');
    case 'TsavgPlot:'
      i_tsplot = i_tsplot + 1;
      Cdata.TsavgPlots(i_tsplot).Var    = Fields{2};
      Cdata.TsavgPlots(i_tsplot).Name   = Fields{3};
      Cdata.TsavgPlots(i_tsplot).Units  = Fields{4};
      Cdata.TsavgPlots(i_tsplot).Title  = Fields{5};
      Cdata.TsavgPlots(i_tsplot).LegLoc = Fields{6};
      Cdata.TsavgPlots(i_tsplot).Ymin   = sscanf(Fields{7},  '%f');
      Cdata.TsavgPlots(i_tsplot).Ymax   = sscanf(Fields{8},  '%f');
    case 'AzavgEofPlot:'
      i_aeofplot = i_aeofplot + 1;
      Cdata.AzavgEofPlots(i_aeofplot).Var   = Fields{2};
      Cdata.AzavgEofPlots(i_aeofplot).Name  = Fields{3};
      Cdata.AzavgEofPlots(i_aeofplot).Units = Fields{4};
      Cdata.AzavgEofPlots(i_aeofplot).Num   = sscanf(Fields{5},  '%d');;
      Cdata.AzavgEofPlots(i_aeofplot).Clim  = sscanf(Fields{6},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Cinc  = sscanf(Fields{7},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Rmin  = sscanf(Fields{8},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Rmax  = sscanf(Fields{9},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Zmin  = sscanf(Fields{10},  '%f');
      Cdata.AzavgEofPlots(i_aeofplot).Zmax  = sscanf(Fields{11}, '%f');
  end
end

end
