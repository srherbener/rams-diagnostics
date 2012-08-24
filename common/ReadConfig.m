function [ Cases, Tdirs ] = ReadConfig ( Cfile )
%ReadConfig read diagnostic configuration file for a set of simulations
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

% for now understand only Case: and TimeDir: keys
% read through entire list of key,data pairs and
%   place the Case: data values into Cases array
%   place the TimeDir: data values into Tdirs array

i_case = 0;
i_tdir = 0;
for i = 1:length(InLines)
  [ Fields ] = Line2Fields(InLines(i,:), ' ');

  if (strcmp(Fields{1},'Case:'))
    i_case = i_case + 1;
    Cases{i_case} = Fields{2};
  else if (strcmp(Fields{1},'TimeDir:'))
    i_tdir = i_tdir + 1;
    Tdirs{i_tdir} = Fields{2};
  end
end

end
