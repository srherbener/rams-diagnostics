%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ] = WriteDpFile(U, V, T, Zg, RH, Lon, Lat, Press, Year, Month, Day, Time, Fid)
% WriteDpFile create a RAMS DP file with initial atmoshperic fields

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


  Nx = length(Lon);
  Ny = length(Lat);
  Nz = length(Press);

  % Header
  fprintf(Fid, '999999 2\n');
  fprintf(Fid, '%d %d %d %d 0 %d %d %d\n', Year, Month, Day, Time, Nz, Nx, Ny);
  fprintf(Fid, '1 %f %f %f %f %f %f 0.0 0.0 0.0\n', Lon(2)-Lon(1), Lat(2)-Lat(1), Lat(1), Lon(1), Lat(end), Lon(end));
  fprintf(Fid, '1');
  for k = 1:Nz
    fprintf(Fid, ' %f', Press(k));
  end
  fprintf(Fid, '\n');

  % Variables
  for k = 1: Nz
    WriteDpVar(U,  Nx, Ny, Nz, k, Fid);
    WriteDpVar(V,  Nx, Ny, Nz, k, Fid);
    WriteDpVar(T,  Nx, Ny, Nz, k, Fid);
    WriteDpVar(Zg, Nx, Ny, Nz, k, Fid);
    WriteDpVar(RH, Nx, Ny, Nz, k, Fid);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ] = WriteDpVar(Var, Nx, Ny, Nz, k, Fid)
% WriteDpVar write out the contents of one level of one varialbe into the DP file

  NumPerLine = 10;

  for j = 1:Ny
    for i = 1:Nx
       fprintf(Fid, ' %f', Var(i,j,k));
    end
    fprintf(Fid, '\n');
  end
end
