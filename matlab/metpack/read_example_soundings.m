function [ P T Td RH Wspd Wdir Nk ] = read_example_soundings(Sfile, Sname, Sset)
%read_example_soundings read in a sounding from examples that came with metpack

% souding input format
% http://weather.uwyo.edu/upperair/sounding.html
%-----------------------------------------------------------------------------
%   PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV
%    hPa     m      C      C      %    g/kg    deg   knot     K      K      K 
%-----------------------------------------------------------------------------

Sstruct = load(Sfile);
S = Sstruct.(Sname).(Sset);

% variables needed by metpack routines
P    = S(:,1);
T    = S(:,3);
Td   = S(:,4);
RH   = S(:,5)/100; % convert percent to decimal fraction
Wdir = S(:,7);
Wspd = S(:,8)/0.514; % convert knots to m/s

Nk = size(S,1); %needed for cape calculation

end
