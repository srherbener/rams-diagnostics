function [ LatKm ] = LatToKm(Lat)
% LatToKm function return latitude in km from latitude given in degrees
%
  R = 6371; % radius of earth in km

  LatRad = Lat * pi / 180;
  LatKm = LatRad * R;
end

