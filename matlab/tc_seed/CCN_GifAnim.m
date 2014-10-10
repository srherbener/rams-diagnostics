% script to create a GIF animation of the contents of the CWP file

clear;

InFileList = { 
  { 'ccn_conc-TCS_GN_C1000-AS-1998-08-22-120000-g3.h5' 'ccn_concen' 0.001 3 55 0 0.1 }
  { 'aerosol_outside_grid3.h5' 'aerosol' 1 3 55 180 0.1 }
  };

OutFileList = {
  'CCN_anim_in_g3.gif'
  'CCN_anim_out_g3.gif'
  };

fprintf('Generating aerosol animations:\n');
fprintf('\n');

Zlev = 1.5; % km

Clevs = (0:20:1000);

for i = 1:2
  InFile  = InFileList{i}{1};
  InVar   = InFileList{i}{2};
  InScale = InFileList{i}{3};
  T1      = InFileList{i}{4};
  T2      = InFileList{i}{5};
  Rotate  = InFileList{i}{6};
  Delay   = InFileList{i}{7};
  
  OutFile = OutFileList{i};
  
  fprintf('****************************************************\n');
  fprintf('Reading: %s (%s)\n', InFile, InVar);
  fprintf('\n');
  
  %CCN_ALLZ = h5read(InFile, InVar);
  %Lon = h5read(InFile, '/x_coords');
  %Lat = h5read(InFile, '/y_coords');
  %Z = h5read(InFile, '/z_coords') * InScale; % convert to km
  
  CCN_DS = ncgeodataset(InFile);
  
  CCN_VAR = CCN_DS.geovariable(InVar);
  LON_VAR = CCN_DS.geovariable('x_coords');
  LAT_VAR = CCN_DS.geovariable('y_coords');
  Z_VAR = CCN_DS.geovariable('z_coords');
  
  % the two files have different levels, grab the closest to Zlev
  Z = Z_VAR.data(:) * InScale;
  Z1 = find(Z <= Zlev, 1, 'last'); % find the last level below Zlev
  
  %CCN is organized as (t,z,y,x);
  CCN = squeeze(CCN_VAR.data(T1:T2,Z1,:,:));
  Nt = size(CCN,1);
  
  %Read in Lon (x) and Lat (y) values
  Lon = LON_VAR.data(:);
  Lat = LAT_VAR.data(:);
  
  C = squeeze(CCN(1,:,:));
  contourf(Lon, Lat, C, Clevs);
  axis tight
  shading flat;
  set(gca,'nextplot','replacechildren','visible','off')
  f = getframe;
  [im,map] = rgb2ind(f.cdata,256,'nodither');
  im(1,1,1,Nt) = 0;
  for k = 1:Nt
    C = squeeze(CCN(k,:,:));
    contourf(Lon, Lat, C, Clevs);
    shading flat;
    f = getframe;
  
    im(:,:,1,k) = imrotate(rgb2ind(f.cdata,map,'nodither'), Rotate);
  end

  fprintf('Writing: %s\n', OutFile);
  fprintf('\n');
  
  imwrite(im,map,OutFile,'DelayTime',Delay,'LoopCount',0);
end




