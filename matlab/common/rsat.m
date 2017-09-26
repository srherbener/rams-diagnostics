function var = rsat( temp, press, vapor, varargin )
% function var = rsat( temp, press, vapor, varargin )
%   written by Leah Grant  09/2013
%
% Calculate saturation mixing ratio w.r.t. liquid in g/kg (default)
% Arguments: temp in K, pressure in hPa, vapor mixing ratio in g/kg
%
% Pass in [] for arguments that aren't needed.
% Default is to calculate saturation mixing ratio w.r.t. liquid; example:
%   rsatw = rsat( temp, press, [] );
%
% Specify a different variable to calculate using argument pairs
%   to varargin.
% Possible pairs (any case is acceptable):
%   'var','rsatw' --> calculate saturation mixing ratio w.r.t. liquid (g/kg)
%                     (this is the default)
%   'var','rsati'  --> calculate saturation mixing ratio w.r.t. ice (g/kg)
%                      temp (K) <= T0
%   'var','esw' --> calculate saturation vapor pressure over liquid (hPa)
%   'var','esi' --> calculate saturation vapor pressure over ice (hPa)
%                   temp (K) <= T0
%   'var','rhw' --> calculate relative humidity w.r.t. liquid (%)
%   'var','rhi' --> calculate relative humidity w.r.t. ice (%)
%                    temp(K) <= T0
%
% Pass in [] for variables that aren't needed.
% For example, if you want to calculate esw, you only need temperature. Thus
%   you can pass in [] for press and vapor. Example call:
%     es = rsat( temp, [], [], 'var', 'esw' );
% esw and esi only need temperature
% rsatw and rsati need temperature and pressure
% rhw and rhi need temperature, pressure, and vapor
%
% Equations used for esw, esi are from Murphy & Koop (2005)


% initialize vartype: 1 for saturation mixing ratio over liquid
vartype = 1;

% get arguments
nArgExtra=length(varargin);
if nArgExtra>0
    for i=1:2:nArgExtra
        if strcmpi(varargin{i},'var')
            varn=lower(varargin{i+1});
            if strcmp(varn,'rsati'); vartype=2;
            elseif strcmp(varn,'rsatw'); vartype=1;
            elseif strcmp(varn,'esw'); vartype=3;
            elseif strcmp(varn,'esi'); vartype=4;
            elseif strcmp(varn,'rhw'); vartype=5;
            elseif strcmp(varn,'rhi'); vartype=6;
            elseif strcmp(varn,'rsat') % other possible inputs: assume liquid
                % but print a message
                disp('Calculating saturation mixing ratio w.r.t. liquid')
                vartype=1;
            elseif strcmp(varn,'es')
                disp('Calculating saturation vapor pressure w.r.t. liquid')
                vartype=3;
            elseif strcmp(varn,'rh')
                disp('Calculating relative humidity w.r.t. liquid')
                vartype=5;
            else
                disp(['Unknown variable type, ',varn,': Default to calculating rsatw'])
            end
        end
    end
end           

T0=273.15;

if any(vartype==[1 3 5]); liquid=1; ice=0;
elseif any(vartype==[2 4 6]); ice=1;
    if any(temp(:)>T0); liquid=1; else liquid=0; end
end

% error checking with arguments
if liquid && ice  % print message
    disp('Variable with respect to ice for temperature > freezing')
    disp('   --> calculated w.r.t. liquid')
end

if any(vartype==[1 2]) && isempty(press)
    error('Need to pass in pressure to calculate sat. mixing ratio!')
end
if any(vartype==[5 6]) && ( isempty(press) || isempty(vapor) )
    error('Need to pass in pressure & vapor mixing ratio to calculate RH!')
end


% define some values
Rd = 287.04; % J/kg/K
Rv = 461.51; % J/kg/K
eps = Rd/Rv;

% compute es for liquid or ice
if liquid
    % Murphy and Koop 2005 eqn. (10); valid for 123<T<332K
    esw = exp( ...
           54.842763 - 6763.22./temp - 4.210*log(temp) + ...
           0.000367*temp + tanh( 0.0415*(temp-218.8) )...
                           .*( 53.878 - 1331.22./temp ...
                               - 9.44523*log(temp) + 0.014025*temp )...
          ); % Pa
end
if ice
    % Murphy and Koop 2005 eqn. (7);  valid for T > 110K
    esi = exp( 9.550426 - 5723.265./temp + ...
              3.53068*log(temp) - 0.00728332*temp ); % Pa
end

if liquid && ~ice
    es = esw;
elseif ~liquid && ice
    es = esi;
elseif liquid && ice
    es = esi;
    es(temp>T0)=esw(temp>T0); % use esw where above freezing
end
es = es/100.; % hPa

% if want to compute esw or esi, we're done..
if any(vartype==[3 4])
    var = es; % hPa
    return % leave function
end

% calculate rsw, rsi 
if any(vartype==[1 2])
    % saturation mixing ratio rsat = eps * es / (p-es)
    var = eps * es ./ ( press - es ); % es, press in hPa
    var = var*1000.; % g/kg
    return
end

% calculate rhw, rhi
if any(vartype==[5 6])
    % rh = e / es
    % calculate e: rv=eps*e/(p-e), so e=rv*p/(eps+rv)
    e = vapor/1000 .* press ./ ( eps+vapor/1000 ); % hPa
    var = e ./ es * 100.; % pct
    return
end