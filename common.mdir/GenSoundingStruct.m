function [ SndStruct ] = GenSoundingStruct(InFilePattern)
% GenSoundingStruct generate a structure with sounding data from the input files given by InFilePattern
%
%   InFilePattern is a string with a file matching pattern that is used to select all intput files
%   SndStruct is a structure that contains all of the soundings found in the input files.
%
%   This function takes all of the input files that match InFilePattern, reads the soundings out of these
%   files and merges them into a structure which is passed back through SndStruct.
%
% Profiles are arranged by height:
%   pressure
%   temp
%   dew point
%   relateve humidity
%   cloud mixing ratio
%   ice mixing ratio
%   rain mixing ratio
%   u wind speed
%   v wind speed

  % Use the pattern to form a list of input files.
  LsInfo = textscan(ls('-1', InFilePattern), '%s');
  InFileList = LsInfo{1};

  % Read in the input files and construct a set of soundings.
  %   Use the heights out of the first file for all sets of profiles
  for i = 1:length(InFileList)
    InFile = InFileList{i};
    fprintf('Reading: %s\n', InFile);

    % The first seven lines are header information:
    %      Line          contents
    %       1             name
    %       2             latitude
    %       3             longitude
    %       4             time step
    %       5             blank
    %       6             column headers
    %       7             blank
    %
    % all remaining lines contain:
    %
    %     10 numbers (%e format): Height, Pressure, Temp, Dew point, RH, Cloud, Ice, Rain, U, V

    fid = fopen(InFile);
    Line = textscan(fid, '%s %s', 1);
    Name = Line{2};

    Line = textscan(fid, '%s %f', 1);
    Lat = Line{2};

    Line = textscan(fid, '%s %f', 1);
    Lon = Line{2};

    Line = textscan(fid, '%s %s %d', 1);
    TimeStep = Line{3};

    Line = textscan(fid, '%s', 20); % grabs the column headers (next 20 tokens, separated by white space)
    Cheaders = Line{1};

    % rest of file is the sounding data
    Line = textscan(fid, '%f %f %f %f %f %f %f %f %f %f');
    if (i == 1)
      SndStruct.Height = [ Line{1} ];
    end
    SndStruct.Pressure{i}    = [ Line{2} ];
    SndStruct.Temperature{i} = [ Line{3} ];
    SndStruct.DewPoint{i}    = [ Line{4} ];
    SndStruct.RelHum{i}      = [ Line{5} ];
    SndStruct.Cloud{i}       = [ Line{6} ];
    SndStruct.Ice{i}         = [ Line{7} ];
    SndStruct.Rain{i}        = [ Line{8} ];
    SndStruct.U{i}           = [ Line{9} ];
    SndStruct.V{i}           = [ Line{10} ];

    fclose(fid);
  end
end
