% script to create a GIF animation of the contents of the CWP file

clear;

% Grab the cloud water path directly from REVU output
CCN_ALLZ = hdf5read('ccn_conc-TCS_GN_C1000-AS-1998-08-22-120000-g3.h5','/ccn_concen');
Lon = hdf5read('ccn_conc-TCS_GN_C1000-AS-1998-08-22-120000-g3.h5','/x_coords')';
Lat = hdf5read('ccn_conc-TCS_GN_C1000-AS-1998-08-22-120000-g3.h5','/y_coords')';

CCN = squeeze(CCN_ALLZ(:,:,17,3:73)); % z = 17 --> ~1500m AGL

DelayTime = 0.1;

Nt = size(CCN,3);

Clevs = (0:20:1000);
C = squeeze(CCN(:,:,1)');
contourf(Lon, Lat, C, Clevs);
axis tight
shading flat;
set(gca,'nextplot','replacechildren','visible','off')
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,Nt) = 0;
for k = 1:Nt
  C = squeeze(CCN(:,:,k)');
  contourf(Lon, Lat, C, Clevs);
  shading flat;
  f = getframe;
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
end
imwrite(im,map,'plots/Asource_movie.gif','DelayTime',DelayTime,'LoopCount',0);
