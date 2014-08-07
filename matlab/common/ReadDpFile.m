%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ U V T Zg RH Lon Lat Press Year Month Day Time ] = ReadDpFile(DpFile)
% ReadDpFile read a RAMS DP file with initial atmoshperic fields, assumes lat/lon projection

  % The file format is described in CSU-RAMS manual:
  %
  %   "CSU-RAMS: Standard Input Formats for Pressure Coordinate and Observational Data (Data Prep Stage)" 
  %
  % written by Stephen M. Saleeby, Feb 12, 2014
  %
  %
  % There are 5 atmospheric variables:
  %    u velocity  (zonal, m/s)
  %    v velocity  (meridonal, m/s)
  %    temperature (K)
  %    geopotential height (m)
  %    relative humidity (fraction)
  %
  % and 6 soil variables:
  %    top soil level volumetric soil moisture (m^3/m^3)
  %    next level down volumetric soil moisture (m^3/m^3)
  %    top soil level soil temperature (K)
  %    next level down soil temperature (K)
  %    snow water equivalent (kg/m^2)
  %    snow depth (m)
  %
  % Okay to skip the soil variables as long as RAMSIN is configured to initialize
  % the soil variables by another method.

  Fid = fopen(DpFile, 'rt');

  % Line 1 - Header, Version
  Header = fscanf(Fid, '%d', 1);  % 999999
  Version = fscanf(Fid, '%d', 1); % 2
  
  % Line 2 - date/time stamp, size info
  Year  = fscanf(Fid, '%d', 1);
  Month = fscanf(Fid, '%d', 1);
  Day   = fscanf(Fid, '%d', 1);
  Time  = fscanf(Fid, '%d', 1);
  Vtime = fscanf(Fid, '%d', 1);
  Nz    = fscanf(Fid, '%d', 1);
  Nx    = fscanf(Fid, '%d', 1);
  Ny    = fscanf(Fid, '%d', 1);

  % Line 3 - description of horizontal coords
  Projection = fscanf(Fid, '%d', 1);
  DeltaLon   = fscanf(Fid, '%f', 1);
  DeltaLat   = fscanf(Fid, '%f', 1);
  SwLat      = fscanf(Fid, '%f', 1);
  SwLon      = fscanf(Fid, '%f', 1);
  NeLat      = fscanf(Fid, '%f', 1);
  NeLon      = fscanf(Fid, '%f', 1);
  TangLat    = fscanf(Fid, '%f', 1);
  CenterLon  = fscanf(Fid, '%f', 1);
  SecondLat  = fscanf(Fid, '%f', 1);

  % reconstruct lat and lon
  Lon = zeros([ Nx 1 ]);
  Lon(1) = SwLon;
  for i = 2:Nx
    Lon(i) = Lon(i-1) + DeltaLon;
  end

  Lat = zeros([ Ny 1 ]);
  Lat(1) = SwLat;
  for i = 2:Ny
    Lat(i) = Lat(i-1) + DeltaLat;
  end

  % Line 4 - vertical coords
  VcoordType = fscanf(Fid, '%d', 1);
  Press = fscanf(Fid, '%f', Nz);

  % Remaining sections are U, V, T, Zg, RH in that order
  % For each variable, data in file is stored in column-major format which
  % matches the order that fscanf will use. However, fscanf can only do 1
  % or 2 dimensional arrays.
  VarSize = [ Nx Ny Nz ];
  U = zeros(VarSize);
  V = zeros(VarSize);
  T = zeros(VarSize);
  Zg = zeros(VarSize);
  RH = zeros(VarSize);

  HorizDims = [ Nx Ny ];
  for k = 1:Nz
    U(:,:,k)  = fscanf(Fid, '%f', HorizDims);
    V(:,:,k)  = fscanf(Fid, '%f', HorizDims);
    T(:,:,k)  = fscanf(Fid, '%f', HorizDims);
    Zg(:,:,k) = fscanf(Fid, '%f', HorizDims);
    RH(:,:,k) = fscanf(Fid, '%f', HorizDims);
  end

  fclose(Fid);

end

%%%   Nx = length(Lon);
%%%   Ny = length(Lat);
%%%   Nz = length(Press);
%%% 
%%%   % Header
%%%   fprintf(Fid, '999999 2\n');
%%%   fprintf(Fid, '%d %d %d %d 0 %d %d %d\n', Year, Month, Day, Time, Nz, Nx, Ny);
%%%   fprintf(Fid, '1 %f %f %f %f %f %f 0.0 0.0 0.0\n', Lon(2)-Lon(1), Lat(2)-Lat(1), Lat(1), Lon(1), Lat(end), Lon(end));
%%%   fprintf(Fid, '1');
%%%   for k = 1:Nz
%%%     fprintf(Fid, ' %f', Press(k));
%%%   end
%%%   fprintf(Fid, '\n');
%%% 
%%%   % Variables
%%%   for k = 1: Nz
%%%     WriteDpVar(U,  Nx, Ny, Nz, k, Fid);
%%%     WriteDpVar(V,  Nx, Ny, Nz, k, Fid);
%%%     WriteDpVar(T,  Nx, Ny, Nz, k, Fid);
%%%     WriteDpVar(Zg, Nx, Ny, Nz, k, Fid);
%%%     WriteDpVar(RH, Nx, Ny, Nz, k, Fid);
%%%   end
%%% 
%%% end
%%% 
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% function [ ] = WriteDpVar(Var, Nx, Ny, Nz, k, Fid)
%%% % WriteDpVar write out the contents of one level of one varialbe into the DP file
%%% 
%%%   NumPerLine = 10;
%%% 
%%%   for j = 1:Ny
%%%     for i = 1:Nx
%%%        fprintf(Fid, ' %f', Var(i,j,k));
%%%     end
%%%     fprintf(Fid, '\n');
%%%   end
%%% end
