%% this is layer interface heights for PW computation
%  always start with this section to initialize
z = [      0.000000000000000000E+00 ...
  0.500000000000000000E+02 ...
  0.101424003601074219E+03 ...
  0.157282501220703125E+03 ...
  0.219504501342773437E+03 ...
  0.288815002441406250E+03 ...
  0.366022003173828125E+03 ...
  0.452024475097656250E+03 ...
  0.547824462890625000E+03 ...
  0.654538940429687500E+03 ...
  0.773410400390625000E+03 ...
  0.905823913574218750E+03 ...
  0.105332287597656250E+04 ...
  0.121762536621093750E+04 ...
  0.140064636230468750E+04 ...
  0.160451782226562500E+04 ...
  0.183161425781250000E+04 ...
  0.208458276367187500E+04 ...
  0.236637060546875000E+04 ...
  0.268026123046875000E+04 ...
  0.302991113281250000E+04 ...
  0.341939453125000000E+04 ...
  0.385324951171875000E+04 ...
  0.433186767578125000E+04 ...
  0.483186767578125000E+04 ...
  0.533186767578125000E+04 ...
  0.583186767578125000E+04 ...
  0.633186767578125000E+04 ...
  0.683186767578125000E+04 ...
  0.733186767578125000E+04 ...
  0.783186767578125000E+04 ...
  0.833186718750000000E+04 ...
  0.883186718750000000E+04 ...
  0.933186718750000000E+04 ...
  0.983186718750000000E+04 ...
  0.103318671875000000E+05 ...
  0.108318671875000000E+05 ...
  0.113318671875000000E+05 ...
  0.118318671875000000E+05 ...
  0.123318671875000000E+05 ...
  0.128318671875000000E+05 ...
  0.133318671875000000E+05 ...
  0.138318671875000000E+05 ...
  0.143318671875000000E+05 ...
  0.148318671875000000E+05 ...
  0.153318671875000000E+05 ...
  0.158318671875000000E+05 ...
  0.163318671875000000E+05 ...
  0.168318671875000000E+05 ...
  0.173318671875000000E+05 ...
  0.178318671875000000E+05 ...
  0.183318671875000000E+05 ...
  0.188318671875000000E+05 ...
  0.193318671875000000E+05 ...
  0.198318671875000000E+05 ...
  0.203359609375000000E+05 ...
  0.208695585937500000E+05 ...
  0.214565156250000000E+05 ...
  0.221021679687500000E+05 ...
  0.228123847656250000E+05 ...
  0.235936230468750000E+05 ...
  0.244529843750000000E+05 ...
  0.253982812500000000E+05 ...
  0.264381093750000000E+05 ...
];

 SimList = {
   { 'RCE_1km'    'RAMS 2mom'      }
   { 'RCE_1km_SM' 'RAMS 1mom'      }
   { 'RCE_1km_DM' 'RAMS 2mom mean' }
   { 'RCE_1km_DP' 'RAMS 2mom base' }
   };
 Nsims = length(SimList);
 
 %user edits fprefix and 
  nx = 480;  % horizontal dimension
  nend=1201; % number of records to process
  dz = z(2:64)-z(1:63);
  n1 = 1;  % record stride

 %% do this section to start a nc file

 for isim = 1:Nsims
   Sim = SimList{isim}{1};
   SimLabel = SimList{isim}{2};

   fprefix = sprintf('RAMS/%s/RAMS/%s-L-2012-', Sim, Sim);
   ncfile = sprintf('%s_pw.nc', Sim);

   fprintf('Simulation: %s\n', Sim);
   fprintf('  File prefix: %s\n', fprefix);
   fprintf('  Output file: %s\n', ncfile);
   fprintf('\n');

   offset = 0;
   ncid = netcdf.create(ncfile,'CLOBBER');
   dimidx = netcdf.defDim(ncid,'x',nx);
   dimidy = netcdf.defDim(ncid,'y',nx);
   NC_UNLIMITED = netcdf.getConstant('NC_UNLIMITED');
   dimidt = netcdf.defDim(ncid,'time',NC_UNLIMITED);
   varidt = netcdf.defVar(ncid,'time','NC_DOUBLE',dimidt);
   varid = netcdf.defVar(ncid,'PW','NC_FLOAT',[dimidx dimidy dimidt]);
   varidlhf = netcdf.defVar(ncid,'LHF','NC_FLOAT',[dimidx dimidy dimidt]);
   varidshf = netcdf.defVar(ncid,'SHF','NC_FLOAT',[dimidx dimidy dimidt]);
   varidpr = netcdf.defVar(ncid,'PREC','NC_FLOAT',[dimidx dimidy dimidt]);
   varido = netcdf.defVar(ncid,'OLR','NC_FLOAT',[dimidx dimidy dimidt]);
   netcdf.endDef(ncid);
    
%     %% do this section to continue an nc file
%     offset = 2352;
%     ncid = netcdf.open(ncfile,'WRITE');
%     varidt = netcdf.inqVarID(ncid,'time');
%     varid = netcdf.inqVarID(ncid,'PW');
%     varidlhf = netcdf.inqVarID(ncid,'LHF');
%     varidshf = netcdf.inqVarID(ncid,'SHF');
%     varidpr = netcdf.inqVarID(ncid,'PREC');
%     varido = netcdf.inqVarID(ncid,'OLR');
     
    %% start time looping the RAMS files and create/extend nc file
    for n = n1+offset:n1:nend;
        
        % adjust from year mon day to record number
        % assume start is 1 Jan 2012
        h = mod(n,24);
        d = 1+(n-h)/24;
        m = 1;
        if n > 743; 
            d = d - 31;
            m = m + 1;
        end;
        if n > 1439; 
            d = d - 29;  !
            m = m + 1;
        end;
        if n > 2183; 
            d = d - 31;
            m = m + 1;
        end;
        if n > 2903; 
            d = d - 30;
            m = m + 1;
        end;
        if n > 3647; 
            d = d - 31;
            m = m + 1;
        end;
        hstring = num2str(h,'%2.2u');
        dstring = num2str(d,'%2.2u');
        mstring = num2str(m,'%2.2u');
        datestring = strcat(mstring,'-',dstring,'-',hstring)
       file=strcat(fprefix,datestring,'0000-g1.h5');
       
       %read in the data
       x = hdf5read(file,'LWUP');
       q = hdf5read(file,'RV');
       rho = hdf5read(file,'DN0');
       sfr = hdf5read(file,'SFLUX_R');
       sft = hdf5read(file,'SFLUX_T');
       p = hdf5read(file,'PCPRR');
       
       %compute precipitable water
       y=zeros(nx,nx);
       q = rho.*q;
       for k=1:1:63;
           y = y + dz(1,k)*q(:,:,k+1);
       end;
       
       %get domain averages, save 2d data in nc file
       pw(n/n1) = mean(mean(y(:,:),2),1);
       evap(n/n1) = mean(mean(sfr(:,:),2),1)*86400;
       fss(n/n1) = mean(mean(sft(:,:),2),1)*1004;
       prec(n/n1) = mean(mean(p(:,:),2),1)*86400;
       olr(n/n1) = double(mean(mean(x(:,:,64),2),1));
       time(n/n1) = n/24;
       netcdf.putVar(ncid,varidt,n/n1-1,1,time(n/n1));
       netcdf.putVar(ncid,varid,[0 0 n/n1-1],[nx nx 1],y);
       netcdf.putVar(ncid,varidlhf,[0 0 n/n1-1],[nx nx 1],sfr*86400);
       netcdf.putVar(ncid,varidshf,[0 0 n/n1-1],[nx nx 1],sft*1004);
       netcdf.putVar(ncid,varidpr,[0 0 n/n1-1],[nx nx 1],p*86400);
       netcdf.putVar(ncid,varido,[0 0 n/n1-1],[nx nx 1],x(:,:,64));
       netcdf.sync(ncid);
    end;
 end

%%% %% Do this section if you are not working with the RAMS h5 files
%%% 
%%% time = ncread(ncfile,'time')';
%%% y = ncread(ncfile,'PW');
%%%  pw(:) = mean(mean(y,2),1);
%%% x = ncread(ncfile,'OLR');
%%%  olr(:) = mean(mean(x,2),1);
%%% p = ncread(ncfile,'PREC');
%%%  prec(:) = mean(mean(p,2),1);
%%% sfr = ncread(ncfile,'LHF');
%%%  evap(:) = mean(mean(sfr,2),1);
%%%  
%%% 
%%% %% time series
%%% sizer = size(time,2);
%%% sizeb=2400; % time series start
%%% sizee = sizer ; % time series end
%%% figure;plot(time(sizeb:sizee),olr(sizeb:sizee));title( strcat('OLR:', SimLabel) );ylabel('W/m^2');xlabel('time, days');
%%% figure;plot(time(sizeb:sizee),pw(sizeb:sizee));title( strcat('Column Water Vapor: ', SimLabel) );ylabel('kg/m^2');xlabel('time, days');
%%% figure;plot(time(sizeb:sizee),prec(sizeb:sizee));title( strcat('Precipitation rate: ', SimLabel) );ylabel('mm/day');xlabel('time, days');
%%% figure;plot(time(sizeb:sizee),evap(sizeb:sizee));title( strcat('SFC Evaporation rate: ', SimLabel) );ylabel('mm/day');xlabel('time, days');
%%% %10 day focus
%%% figure;plot(time(sizee-239:sizee),prec(sizee-239:sizee));title( strcat('Precipitation rate: ', SimLabel) );ylabel('mm/day');xlabel('time, days');
%%% 
%%% %% power spectra
%%% sizer = size(time,2);
%%% sizee = sizer ; % time series end
%%% % work with 60 days
%%% powerspectraplot(prec(sizee-1440:sizee));title(strcat('Precipitation spectra: ', SimLabel ));
%%% powerspectraplot(pw(sizee-1440:sizee));title(strcat('Precipitable water spectra: ', SimLabel ));
%%% powerspectraplot(olr(sizee-1440:sizee));title(strcat('OLR spectra: ', SimLabel ));
%%% powerspectraplot(evap(sizee-1440:sizee));title(strcat('Evaporation spectra: ', SimLabel ));
