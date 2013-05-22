function [ Profile ] = SoundingToProfile(InFilePattern)
% SoundingToProfile generate a profile structure from a list of files created by DumpSounding.gs
%
%   InFilePattern is a string with wildcards that is used to select intput files for the sounding
%   Profile is a structure with selected pieces from Sounding
%
% For now, the following are placed in Profile
%
%   Height
%   Temperature
%   Relative Humidity
%   Cloud water
%   Cloud ice
%   Precip (rain)
%

  Snames = {
    'Height'
    'Temperature'
    'RelHum'
    'Cloud'
    'Ice'
    'Rain'
    };

  Pnames = {
    'Heights'
    'Temperatures'
    'RelHum'
    'CloudWater'
    'CloudIce'
    'Precip'
    };

  % Collect the sounding data and copy over the desired structures
  S = GenSoundingStruct(InFilePattern);

  % Use the mapping given in Snames and Pnames to do the copy
  for i = 1:length(Snames)
    Profile.(Pnames{i}) = S.(Snames{i});
  end

end
