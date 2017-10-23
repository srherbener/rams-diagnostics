function [ LonKm ] = LonToKm(Lon, Lat)
% LonToKm function return longitude in km from longitude in degrees, and reference latitude in degrees
%
  R = 6371; % radius of earth in km

  LonRad = Lon * pi / 180;
  LatRad = Lat * pi / 180;

  LonKm = LonRad * R * cos(LatRad);
end

