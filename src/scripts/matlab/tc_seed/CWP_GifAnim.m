% script to create a GIF animation of the contents of the CWP file

clear;

% Grab the cloud water path directly from REVU output
CWP = hdf5read('vint_cloud-TCS_CNTL-AS-1998-08-22-120000-g3.h5','/vertint_cloud');
Lon = hdf5read('vint_cloud-TCS_CNTL-AS-1998-08-22-120000-g3.h5','/x_coords')';
Lat = hdf5read('vint_cloud-TCS_CNTL-AS-1998-08-22-120000-g3.h5','/y_coords')';

Nt = size(CWP,3);

Clevs = (0:0.2:2);
C = squeeze(CWP(:,:,1)');
contourf(Lon, Lat, C, Clevs);
axis tight
shading flat;
set(gca,'nextplot','replacechildren','visible','off')
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,Nt) = 0;
for k = 1:Nt
  C = squeeze(CWP(:,:,k)');
  contourf(Lon, Lat, C, Clevs);
  shading flat;
  f = getframe;
  im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
end
imwrite(im,map,'CWP_anim.gif','DelayTime',0,'LoopCount',0);