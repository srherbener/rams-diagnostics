function zlevs = read_zlevs_hfile( directory, levels )
% function zlevs = read_zlevs_hfile( dir, levels )
% 
% reads z-levels from the first header file in the supplied directory
% 
% output: zlevs in m
%
% input: 
%   dir (required) - directory where header files are listed
%   levels (optional) - either "t" for scalar levels or "m" for momentum levels
%     default is to read t levels
%

% check arguments
if nargin==0
    error('Please input the directory name containing header files as the first argument')
elseif nargin<=2 
    % check directory argument
    if ~isdir(directory)
        error(['Directory ',directory,' does not exist'])
    elseif isempty(dir([directory,'/*head.txt']))
        error(['No header files in directory ',directory])
    end
    if nargin==1
        levels='t'; % set second argument to default
    elseif ~any(strcmp(levels,{'t';'m'}))
        error('Second argument "levels" should be either "t" or "m"')
    end
else
    error('Too many input arguments')
end
    
% set z levels string to search for in header file
if strcmp(levels,'t'); zread='ztn01'; disp('Reading t levels')
elseif strcmp(levels,'m'); zread='zmn01'; disp('Reading m levels')
end

files=dir([directory,'/*head.txt']);
hfile=[directory,files(1).name]; clear files % get the first header file

% get the # of header lines by finding the line where z string is printed
[~,nhlines]=system(['grep -m 1 -n ',zread,' ',hfile]);
i=find(nhlines==':'); % only use first part of grep output: line number
nhlines=str2num(nhlines(1:i-1))+1; % add 1 b/c line below the z variable string 
 % is the number of levels, so we want to skip this line

zlevs = importdata(hfile,' ',nhlines); zlevs=zlevs.data; % z levels array in m