function [ Exps, Tdirs ] = ReadConfig ( Cfile )
%ReadConfig read diagnostic configuration file for a set of simulations
%
% This function will read a simple configuration file for processing
% a set of experimental simulation runs.
%
% The config file has the format:
%    Exp: <directory>
%    Tdir: <directory>
%
% Exp: and Tdir: define all of the directories where hdf5 files live
% 
% Exp: defines separate experiments (simulations)
%
% Tdir: defines multiple directories under a single experiment where
%       hdf5 files live
%
% If you have
%   Exp: RUN1
%   Exp: RUN2
%   Exp: RUN3
%   Tdir: T1
%   Tdir: T2
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

% read in the key, data pairs in the config file

[ ckey, cdata ] = textread(char(Cfile), '%s %s');

% for now understand only Exp: and Tdir: keys
% read through entire list of key,data pairs and
%   place the Exp: data values into Exps array
%   place the Tdir: data values into Tdirs array

i_exp = 0;
i_tdir = 0;
for i = 1:length(ckey)
  if (strcmp(ckey(i),'Exp:'))
    i_exp = i_exp + 1;
    Exps(i_exp) = cdata(i);
  else if (strcmp(ckey(i),'Tdir:'))
    i_tdir = i_tdir + 1;
    Tdirs(i_tdir) = cdata(i);
  end
end

end
