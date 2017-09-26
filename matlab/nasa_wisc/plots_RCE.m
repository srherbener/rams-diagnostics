function [] = plots_RCE(ttlname, ftype, iread)


% contents of Leah's startup.m file on frost
% set default line color order for line plots using "hold all"
DefaultColorOrder = ...
[     0.00      0.00      0.00  % black
    1.0000         0         0  % red
         0    0.7500    0.7500  % cyan
         0         0    1.0000  % blue
         0    0.5000         0  % darker green
    0.7500    0.7500         0  % gold
    0.7500         0    0.7500];% pinkish

% default line color and style order
set(0,'DefaultAxesColorOrder',DefaultColorOrder)
set(0,'DefaultAxesLineStyleOrder','-|--|-.|:')

% default axes font size and line width
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultLineLineWidth',4)




%% edit these!
% Note: changed 2D rad top level vars to the 3D ones...
%clear
% directory with RAMS output files
%ttlname='channel-dp-u7-1hd18';
%ttlname='RCE_3km_1mom';
ttlname2=ttlname; ttlname2(ttlname=='_')='-';
%ftype='L'; % "A" for analysis or "L" for lite files
%iread=1; % read files? 1=read, 0=load matlab file

%RAMSdir = ['/Volumes/avalanche/ldgrant/RCE_YS/z.',ttlname,'/NOBAK/'];
%RAMSdir = ['/Volumes/avalanche/ldgrant/rams_RCE_20150928_dev/test/z.RCE-',ttlname,'/NOBAK/'];
RAMSdir = ['/Volumes/avalanche/sherbener/NOBAK/projects/NasaWisc/',ttlname,'/RAMS/'];
%SaveDir = ['./GROW/RCE-',ttlname,'/',ftype,'files/'];
SaveDir = ['./RCE/',ttlname,'/',ftype,'files/'];

SaveGoogleDrive=0; % copy .png files and vintcond directory to google drive?
GoogleDrivePathBase='~/Google Drive/work/RCE plots/';
GoogleDrivePath=[GoogleDrivePathBase,ttlname,'/'];

% grid setup
dx = 3.; % horiz grid spacing in km
%nx = 64; ny = 64; % # horizontal grid points
nx = 256; ny = 256; % # horizontal grid points

% grid point indices for reading and plotting
% x=1 and x=2 are the same; x=end-1 and x=end are the same, and same
% as x=1 and 2 -- so, 3:end-1 are unique x's; same for y I think (?)
x1=3; x2=nx-1;
y1=3; y2=ny-1;
z1=2; z2=150;
tskip = 1; % skip interval to read
tend = 9999; % last file to read


%% set up some arrays/values based on above settings
x = [0:nx-1]*dx; y = [0:ny-1]*dx; % km
zt = read_zlevs_hfile( RAMSdir, 't' )/1000.; % km
zm = read_zlevs_hfile( RAMSdir, 'm' )/1000.;
nz = length(zt);
if z2>nz; z2=nz; end
if ny==1; y1=1; y2=1; end % for 2D

% spatial arrays for reading and plotting (km)
xp = x(x1:x2); nxp=length(xp);
yp = y(y1:y2); nyp=length(yp);
zp = zt(z1:z2); nzp=length(zp);

% files
files = dir([RAMSdir,'/R*-',ftype,'-*h5']); % TEMP added R out front

% time info
nt = min(tend,length(files)); % # output times
%for t=1:tskip:nt; timestr(t,:)=files(t).name(5:21); end % yyyy-mm-dd-HHMMSS
for t=1:tskip:nt; timestr(t,:)=files(t).name(findstr(files(t).name,['-',ftype,'-'])+[3:19]); end % yyyy-mm-dd-HHMMSS
timevec = datevec(timestr,'yyyy-mm-dd-HHMMSS'); % Date Vector array
time1 = repmat(timevec(1,:),[size(timevec,1) 1]); 
tp = etime(timevec,time1) / 3600/24.; % time plotting array (days)
ntp=length(tp);
% dt represents time in days between current file and the file prior; dt(1)=dt(2)
dt = zeros(size(tp)); dt(2:end)=diff(tp); dt(1)=dt(2);

% make output dir if it doesn't exist
if exist(SaveDir,'dir')~=7; mkdir(SaveDir); end

% save dims
save([SaveDir,'dims'],'xp','yp','zp','tp')

%% theta - domain-average theta in time; also do rv and RH
if iread;
    n=0;
    for t=1:tskip:nt; n=n+1;
        theta = h5read( [RAMSdir,files(t).name],'/THETA', [x1 y1 z1], [nxp nyp nzp] );
        exner = h5read( [RAMSdir,files(t).name], '/PI', [x1 y1 z1], [nxp nyp nzp] ); % full exner * cp
        meantheta(:,n) = mean(mean(theta,2),1);
        theta = theta.*exner/1004.; % now this is temp(K)
        meantemp(:,n) = mean(mean(theta,2),1);
        rv = h5read( [RAMSdir,files(t).name], '/RV', [x1 y1 z1], [nxp nyp nzp] )*1000.; % g/kg
        meanrv(:,n) = mean(mean(rv,2),1);
        exner = 1000.*(exner/1004.).^(1004./287.); % now this is press (mb)
        RH = rsat( theta, exner, rv, 'var','rhw' ); % RH (%)
        meanRH(:,n) = mean(mean(RH,2),1);
        if mod(tp(t),5)==0; fprintf('\n'); disp(['day ',num2str(tp(t))]); end
        fprintf('.')
    end
    save([SaveDir,'meanTs'],'meantheta','meantemp')
    save([SaveDir,'meanvapors'],'meanrv','meanRH')
    clear theta exner rv RH
else
    load([SaveDir,'meanTs']); load([SaveDir,'meanvapors']);
end

for ivar=1:3
    figure; set_figsize(gcf,[15 8])
    if ivar==1; pvar=meantheta; varclim=[-15 15]; varttl='\theta'; varunit='K'; varprt='theta';
    elseif ivar==2; pvar=meanrv; varclim=[-4 4]; varttl='r_v'; varunit='g/kg'; varprt='rv';
    elseif ivar==3; pvar=meanRH; varclim=[-50 50]; varttl='RH'; varunit='%'; varprt='RH';
    end
    contourf(tp,zp,bsxfun(@minus,pvar,pvar(:,1)),32,'linestyle','none'); colorbar
    set(gca,'clim',varclim)
    title(['Domain-mean ',varttl,' difference from day 0 (',varunit,') - ',ttlname2])
    %load_cmaps; set(gcf,'colormap',redblu)
    xlabel('time (days)'); ylabel('height (km)')
    print('-dpng',[SaveDir,'avg',varprt,'_hov'])
    
    if ivar==1; pvar=meantemp; varttl='temp'; varprt='temp';
        var0plot=250; varfact1plot=10;
    elseif ivar==2; var0plot=0; varfact1plot=1;
    elseif ivar==3; var0plot=0; varfact1plot=1;
    end
    figure; set_figsize(gcf,[10 10])
    plot(pvar(:,1),zp); hold all
    plot(pvar(:,end),zp);
    plot((pvar(:,end)-pvar(:,1))*varfact1plot+var0plot,zp)
    if ivar==2; ylim([zp(1) 15]); else ylim(zp([1 end])); end
    ylabel('height (km)')
    title(['Domain-mean ',varttl,' (',varunit,') - ',ttlname2])
    if ivar==1; lgd3=['difference*',num2str(varfact1plot),' +',num2str(var0plot)];
    else lgd3='difference';
    end
    legend('day 0',['day ',num2str(tp(end))],lgd3,'location','northeast')
    plot([0 0]+var0plot,zp([1 end]),'k-','linewidth',2)
    if ivar==1
        plot([0 0]+var0plot+varfact1plot,zp([1 end]),'k-','linewidth',1)
        plot([0 0]+var0plot-varfact1plot,zp([1 end]),'k-','linewidth',1)
    end
    grid on
    disp(['Writing: ',SaveDir,'avg',varprt,'_prof']);
    print('-dpng',[SaveDir,'avg',varprt,'_prof'])
end

%% domain-average theta at one time difference from base state
% also average theta over clear regions difference from base state
% plots for Sue

if 0
    % base state from time 0; just read one column
    theta0 = squeeze( h5read( [RAMSdir,files(1).name],'/THETA', [x1 y1 z1], [1 1 nzp] ) );
    
    figure; hold all; startup % add plots as loop goes on
    il=0; % initialize legend counter
    PWthresh = 55; % mm
    ils=0; lstyles{1}='-'; lstyles{2}=':';
    for t=[180 355];
        theta = h5read( [RAMSdir,files(t).name],'/THETA', [x1 y1 z1], [nxp nyp nzp] );
        meantheta = squeeze( mean(mean(theta,1),2) );
        % pull vert int vapor from REVU file
        %PW = h5read( [REVUdir,'vint_vapor-RCE_S300-AC-2012-01-01-000000-g1.h5'],'/vertint_vapor', [x1 y1 t], [nxp nyp 1] );
        % read PW from saved file
        load([SaveDir,'planviews.mat']);
        % get mean theta where pw<thresh
        [meanthetadry,meanthetamoist] = deal( zeros(size(meantheta)) );
        for z=1:nzp
            thtemp = theta(:,:,z);
            meanthetadry(z) = mean(mean(thtemp(PW(:,:,t)<PWthresh),1),2);
            meanthetamoist(z) = mean(mean(thtemp(PW(:,:,t)>=PWthresh),1),2);
        end
        clear thtemp
        
        il=il+1; lgd{il}=['day ',num2str(tp(t)),', domain-avg'];
        il=il+1; lgd{il}=['day ',num2str(tp(t)),', avg PW<',num2str(PWthresh),'mm'];
        il=il+1; lgd{il}=['day ',num2str(tp(t)),', avg PW>',num2str(PWthresh),'mm'];
        
        ils=ils+1; % line style counter
        plot(meantheta-theta0,zp,lstyles{ils},'color',DefaultColorOrder(1,:))
        plot(meanthetadry-theta0,zp,lstyles{ils},'color',DefaultColorOrder(2,:))
        plot(meanthetamoist-theta0,zp,lstyles{ils},'color',DefaultColorOrder(3,:))
    end
    legend(lgd,'location','best')
    ylim([0 25]); ylabel('height (km)')
    title(['\theta differences from base state (K) - ',ttlname2])
    plot([0 0],get(gca,'ylim'),'k-','linewidth',1)
end


%%  equilibrium moisture calculations
if iread;
    n=0;
    for t=1:tskip:nt; n=n+1;
        %accpr(:,:,n) = h5read( [RAMSdir,files(t).name],'/ACCPR', [x1 y1], [nxp nyp] );% ...
            %+ h5read( [RAMSdir,files(t).name],'/ACCPH', [x1 y1], [nxp nyp] ) ...
            %+ h5read( [RAMSdir,files(t).name],'/ACCPD', [x1 y1], [nxp nyp] ); % kg/m2
        LHF = h5read( [RAMSdir,files(t).name], '/SFLUX_R', [x1 y1], [nxp nyp] )*2.5e6; % W/m2
        meanLHF(n) = mean(LHF(:));
        % TEMP: since Steve's files don't have ACCPR
        pcp = h5read( [RAMSdir,files(t).name], '/PCPRR', [x1 y1], [nxp nyp] )*3600*24.; % kg/m2/day
        meanpcp(n) = mean(pcp(:));
        %if n==1
        %    meanpcp(n) = 0;
        %else
        %    % difference = ~average pcp rate over X days
        %    meanpcp(n) = mean(mean(accpr(:,:,n)-accpr(:,:,n-1),2),1)/dt(n);
        %end
        if mod(tp(t),5)==0; disp(['day ',num2str(tp(t))]); end
    end
    save([SaveDir,'sfc-moisture'],'meanpcp','meanLHF')
else
    load([SaveDir,'sfc-moisture']);
end

startup; cols=DefaultColorOrder;
% 5-day moving averages for the smallest output interval; this means that
% if output interval changes, the 5-day moving average will actually be
% more than a 5-day moving average where dt is coarser
daymean = 5;
avgdt = ceil( daymean/min(dt) );
avgdt = avgdt + mod(avgdt+1,2); % this makes sure avgdt is an odd number

figure; set_figsize(gcf,[15 8]);
plot(tp,meanLHF,'color',cols(1,:)); hold all
plot(tp,meanpcp*2.5e6/(24*3600),'color',cols(3,:)) % kg/m2/day*J/kg*days/s = J/m2/s=W/m2
plot(tp,smooth121(meanLHF,avgdt),'linewidth',2,'color',cols(6,:))
plot(tp,smooth121(meanpcp*2.5e6/(24*3600),avgdt),'linewidth',2,'color',cols(4,:))
xlim(tp([1 end])); ylim([0 200]); grid on
xlabel('time (days)')
title(['Domain-mean moisture quantities (W/m^2) - ',ttlname2])
legend('LHF','Precip',['LHF ',num2str(daymean),'-day mean'],...
    ['Precip ',num2str(daymean),'-day mean'],'location','southeast')
disp(['Writing: ', SaveDir,'equil_moisture']);
print('-dpng',[SaveDir,'equil_moisture'])

% animation of precip
if 0
    if tskip==1
        if exist([SaveDir,'accpcp'],'dir')~=7; mkdir([SaveDir,'accpcp']); end
        figure; set_figsize(gcf,[15 5]); %[10 10]);
        for t=2:nt
            clf
            pcolor(xp,yp,(accpr(:,:,t)-accpr(:,:,t-1))'/dt(t)); shading flat; colorbar; set(gca,'clim',[0 150])
            xlabel('x (km)'); ylabel('y (km)')
            title(['Accumulated Precip over ',num2str(dt(t)),sprintf(' days (mm) - day %6.2f',tp(t))])
            print('-dpng',[SaveDir,'accpcp/',sprintf('accpcp_day%06.2f',tp(t)),'.png'])
            pause(.1)
        end
    else
        disp('set tskip to 1 to make precip animation')
    end
end


%% equilibrium energy calculations
% first try to make plots that Herbie made
% plot total surface energy flux (SHF + LHF) and Net Rad Flux Divergence
% for equilibrium: should have net rad cooling of atm ~ sfc fluxes (energy out ~ energy in);

if iread;
    n=0;
    for t=1:tskip:nt; n=n+1;
        fn = [RAMSdir,files(t).name]; % file name
        fn0 = [RAMSdir,files(1).name]; % file name
        
        % sfc fluxes
        SHF = h5read( fn, '/SFLUX_T', [x1 y1], [nxp nyp] )*1004.; % W/m2
        LHF = h5read( fn, '/SFLUX_R', [x1 y1], [nxp nyp] )*2.5e6; % W/m2
        
        % radiation variables
        % The swdn variable at very top level is the same as at top-1 and
        %  as rswdtop; same with swup at top, top-1, and rswutop, and
        %  lwup at top, top-1, and rlontop
        % This is because in rad_mclat.f90, extra rad levels are only used when top press
        %  level > 30mb (3000 Pa); in RCE sims, p at level 75 (25.2km) is 2348Pa
        % NOTE: 3D rad vars are on m levels
        SWdnSfc = h5read( fn, '/RSHORT', [x1 y1], [nxp nyp] );% sfc downwelling SW
        LWdnSfc = h5read( fn, '/RLONG', [x1 y1], [nxp nyp] ); % sfc downwelling LW
        %LWupSfc = h5read( fn, '/RLONGUP', [x1 y1], [nxp nyp] ); % sfc upwelling LW
        LWupSfc = h5read( fn, '/LWUP', [x1 y1 1], [nxp nyp 1] ); % sfc upwelling LW - not correct from the 2D var
        Albedo = h5read( fn, '/ALBEDT', [x1 y1], [nxp nyp] ); % SWupSfc = a*SWdnSfc
        %SWdnTop = h5read( fn, '/RSWDTOP', [x1 y1], [nxp nyp] ); % these 3 vars from Sue
        %SWupTop = h5read( fn, '/RSWUTOP', [x1 y1], [nxp nyp] );
        %LWupTop = h5read( fn, '/RLONTOP', [x1 y1], [nxp nyp] );
        SWdnTop = h5read( fn, '/SWDN', [x1 y1 nz-1], [nxp nyp 1] );
        SWupTop = h5read( fn, '/SWUP', [x1 y1 nz-1], [nxp nyp 1] );
        LWupTop = h5read( fn, '/LWUP', [x1 y1 nz-1], [nxp nyp 1] );
        % net rad cooling = net out (upward) at top + net out (downward) at sfc (W/m2)
        NetRad = +( LWupTop + SWupTop - SWdnTop )... % net rad out at top; LWdnTop=0
            +( -LWupSfc + LWdnSfc + (1-Albedo).*SWdnSfc ); % net rad out at sfc
        % at sfc, swdn-swup = swdn-a*swdn = (1-a)*swdn
        
        % average rad cooling and sfc fluxes
        meanNetRad(n) = mean(NetRad(:));
        meanSFlux(n) = mean(SHF(:)+LHF(:));
        
        % save individual components
        meanSHF(n) = mean(SHF(:));
        meanLHF(n) = mean(LHF(:));
        meanSWdnSfc(n) = mean(SWdnSfc(:));
        meanSWupSfc(n) = mean(Albedo(:).*SWdnSfc(:));
        meanLWdnSfc(n) = mean(LWdnSfc(:));
        meanLWupSfc(n) = mean(LWupSfc(:));
        meanSWdnTop(n) = mean(SWdnTop(:));
        meanSWupTop(n) = mean(SWupTop(:));
        meanLWupTop(n) = mean(LWupTop(:));
        
        % net integrated radiative flux divergence from rad cooling profiles
        % from Manabe and Strickler:  cp/g * int[Psfc->Ptop](dTrad/dt dp) = J/m2/s
        % OR: cp * int[0->ztop](rho*dTrad/dt dz)  checked: these give the same results
        cp=1004; Rd=287;
        % radiative heating rate on t levels. Note: this is a theta-il
        % tendency. So treat as a theta tendency (tried correction for theta-il
        % to theta but that made < 1% difference and anyway isn't in the code)
        RadHtRate = h5read( fn, '/FTHRD', [x1 y1 2], [nxp nyp nz-2] ); % read from 2:nz-1
        % modify RadHtRate so instead of d(thetail)/dt it is d(temp)/dt
        exner = h5read( fn, '/PI', [x1 y1 2], [nxp nyp nz-2] )/cp; % non-dimensional exner
        RadHtRate = RadHtRate.*exner;
        % get base state density
        exner0 = h5read( fn0,'/PI', [x1 y1 2], [nxp nyp nz-2] )/cp;
        press0 = 100000.*exner0.^(cp/Rd); % pressure (Pa)
        Tvirt0 = ( exner0 .* h5read( fn0,'/THETA',[x1 y1 2],[nxp nyp nz-2] ) ) ...
            .* ( 1.+0.61*h5read( fn0,'/RV',[x1 y1 2],[nxp nyp nz-2] ) );
        rho0 = press0/Rd./Tvirt0; clear press0 Tvirt0 exner0
        % get 3D height grid
        [~,~,zzm]=meshgrid(yp,xp,zm(1:nz-1)); % zm=1:nz-1 gives rad cooling rate at zt=2:nz-1
        % integrate: [zm2-zm1] * RadHtRate(zt2) * rho0(zt2) + ...
        IntRadCool = cp * sum( -RadHtRate .* rho0 .* diff(zzm,1,3)*1000, 3 );
        
%         % check 3D rad cooling rates
%         RadCool = -cp*RadHtRate.*rho0.*diff(zzm,1,3)*1000; % W/m2
%         RadCool = -RadHtRate; % K/s
%         
%         lwup=h5read(fn,'/LWUP',[x1 y1 1],[nxp nyp nz-1]);
%         lwdn=h5read(fn,'/LWDN',[x1 y1 1],[nxp nyp nz-1]);
%         swup=h5read(fn,'/SWUP',[x1 y1 1],[nxp nyp nz-1]);
%         swdn=h5read(fn,'/SWDN',[x1 y1 1],[nxp nyp nz-1]);
%         % sfc values aren't right
%         lwup(:,:,1)=h5read(fn,'/RLONGUP',[x1 y1],[nxp nyp]);
%         lwdn(:,:,1)=h5read(fn,'/RLONG',[x1 y1],[nxp nyp]);
%         swdn(:,:,1)=h5read(fn,'/RSHORT',[x1 y1],[nxp nyp]);
%         swup(:,:,1)=swdn(:,:,1).*h5read(fn,'/ALBEDT',[x1 y1],[nxp nyp]);
%         RadCool2 = -diff(swdn,1,3)+diff(swup,1,3)-diff(lwdn,1,3)+diff(lwup,1,3);
%         RadCool2 = RadCool2./(cp*rho0.*diff(zzm,1,3)*1000);
%         % RadCool and RadCool2 match very closely to within ~10^-5 W/m2
%         % RLONGUP and LWUP at z=1 still do not match - 
%         % lwup at z=1 is ~5.2 W/m2 less than RLONGUP - constant across the domain
%         % this is because subroutine sfcrad, which is called in radiation
%         % and in LEAF, is called every timestep for LEAF - so rlongup,
%         % even if overwritten to be consistent with what radiation scheme
%         % computes for lwup(z=1), is overwritten again from LEAF call.
%         % *** This ~5 W/m2 lower for lwup(1) from rad scheme amounts to almost 1K cooler SST ***
%         % swdn at z=1 is identical to RSHORT
%         % lwdn at z=1 is identical to RLONG
%         % swup at z=1 is within 10^-6 of ALBEDT*RSHORT
%         % swdn at nz and nz-1 is identical to RSWDTOP
%         % swup at nz and nz-1 is identical to RSWUTOP
%         % lwup at nz and nz-1 is identical to RLONTOP
        
        meanIntRadCool(n) = mean(IntRadCool(:));
        
        if mod(tp(t),5)==0; disp(['day ',num2str(tp(t))]); end
    end
    save([SaveDir,'EnergyBalance'],'meanSFlux','meanNetRad','meanIntRadCool',...
        'meanSHF', 'meanLHF', 'meanSWdnSfc', 'meanSWupSfc', 'meanLWdnSfc',...
        'meanLWupSfc', 'meanSWdnTop', 'meanSWupTop', 'meanLWupTop')
else
    load([SaveDir,'EnergyBalance']);
end

% 5-day moving averages for the smallest output interval; this means that
% if output interval changes, the 5-day moving average will actually be
% more than a 5-day moving average where dt is coarser
daymean = 5;
avgdt = ceil( daymean/min(dt) );
avgdt = avgdt + mod(avgdt+1,2); % this makes sure avgdt is an odd number
startup; cols=DefaultColorOrder;

% mean net energy balances
figure; set_figsize(gcf,[15 8]);
plot(tp,meanSFlux,'color',cols(2,:)); hold all;
%plot(tp,meanIntRadCool,'color',cols(1,:))
plot(tp,meanNetRad,'color',cols(3,:))
plot(tp,smooth121(meanSFlux,avgdt),'color',cols(5,:),'linewidth',2)
%plot(tp,smooth121(meanIntRadCool,avgdt),'color',cols(6,:),'linewidth',2)
plot(tp,smooth121(meanNetRad,avgdt),'color',cols(4,:),'linewidth',2)
xlim(tp([1 end])); ylim([-20 150]); grid on
xlabel('time (days)'); title(['Domain-mean energy fluxes (W/m^2) - ',ttlname2])
%legend('Sfc Fluxes','Net Int Rad Cool','Net Rad Flux Div','Sfc Fluxes 5-day mean',...
%    'Net Int Rad Cool 5-day mean','Net Rad Flux Div 5-day mean','location','northeast')
legend('Sfc Fluxes','Net Rad Cooling','Sfc Fluxes (running mean)',...
    'Net Rad Cooling (running mean)','location','southeast')
disp(['Writing: ', SaveDir,'equil_energy_net']);
print('-dpng',[SaveDir,'equil_energy_net'])

% mean net sfc flux components
figure; set_figsize(gcf,[15 8]);
plot(tp,meanLHF); hold all;
plot(tp,meanSHF);
xlim(tp([1 end])); ylim([0 100]); grid on
xlabel('time (days)'); title(['Domain-mean sfc fluxes (W/m^2) - ',ttlname2])
legend('LHF','SHF','location','best')
disp(['Writing: ', SaveDir,'equil_energy_sfcFluxes']);
print('-dpng',[SaveDir,'equil_energy_sfcFluxes'])

% mean rad components
figure; set_figsize(gcf,[15 8]);
plot(tp,meanSWdnSfc,'color',cols(7,:)); hold all;
plot(tp,meanSWupSfc,'color',cols(2,:));
plot(tp,meanLWdnSfc,'color',cols(3,:));
plot(tp,meanLWupSfc,'color',cols(1,:));
plot(tp,meanSWdnTop,':','color',DefaultColorOrder(7,:))
plot(tp,meanSWupTop,':','color',DefaultColorOrder(2,:))
plot(tp,meanLWupTop,':','color',DefaultColorOrder(1,:))
xlim(tp([1 end])); ylim([0 500]); grid on
xlabel('time (days)'); title(['Domain-mean rad fluxes (W/m^2) - ',ttlname2])
legend('SWdn sfc','SWup sfc','LWdn sfc','LWup sfc','SWdn top','SWup top','LWup top','location','best')
disp(['Writing: ', SaveDir,'equil_energy_rad']);
print('-dpng',[SaveDir,'equil_energy_rad'])


%% domain-average aerosol concentration

if 0
if iread;
    n=0;
    for t=1:tskip:nt; n=n+1;
        cccnp = h5read( [RAMSdir,files(t).name],'/CCCNP', [x1 y1 z1], [nxp nyp nzp] );
        meanccn(:,n) = mean(mean(cccnp,2),1);
        if mod(tp(t),5)==0; disp(['day ',num2str(tp(t))]); end
    end
    save([SaveDir,'meanccn'],'meanccn')
    clear cccnp
else
    load([SaveDir,'meanccn']);
end

figure; set_figsize(gcf,[15 8]);
contourf(tp,zp,meanccn*1.e-6,32,'linestyle','none'); colorbar
xlabel('time (days)')
ylabel('height (km)')
title(['Domain-mean CCN (#/mg) - ',ttlname2])

figure; set_figsize(gcf,[10 10]);
plot(meanccn(:,1)*1.e-6,zp); hold all
plot(meanccn(:,end)*1.e-6,zp);
set(gca,'ylim',zp([1 end]))
ylabel('height (km)')
title(['Domain-mean CCN (#/mg) - ',ttlname2])
legend('day 0',['day ',num2str(tp(end))],'location','best')
grid on
end

%% check some cross sections and profiles
if 0
    
    t=nt; xplot=50; yplot=50;
    theta = h5read( [RAMSdir,files(t).name],'/THETA' );
    rv = h5read( [RAMSdir,files(t).name],'/RV' )*1000.; % g/kg
    [~,~,N2] = gradient(theta.*(1+0.61*rv/1000.),y,x,zt*1000); % use theta-v
    N2 = 9.8./(theta.*(1+0.61*rv/1000.)).*N2; % g/thetav*dthetav/dz - moist Brust-Vaisala freq.
    
    totcond = ( h5read( [RAMSdir,files(t).name],'/RTP' ) - h5read( [RAMSdir,files(t).name],'/RV' ) )*1000.;
    
    w = h5read( [RAMSdir,files(t).name],'/WP' );
    u = h5read( [RAMSdir,files(t).name],'/UP' );
    press = h5read( [RAMSdir,files(t).name],'/PI' )/1004; % dimensionless
    tempK = theta.*press; % K [here press is still pi]
    press = 1000.*press.^(1004/287); % mb 
    RH = rsat( tempK, press, rv, 'var','rhw' ); % pct
    clear press tempK
    
    % cross sections
    for n=1:6
        figure;
        if n==1; contourf(x,zt,squeeze(N2(:,yplot,:))',32,'linestyle','none'); set(gca,'clim',[-1 10]*1e-4)
            title(['N^2 (s^{-2}), day ',num2str(tp(t)),', through y=',num2str(yplot)])
        elseif n==2; contourf(x,zt,squeeze(totcond(:,yplot,:))',32,'linestyle','none'); set(gca,'clim',[0 5])
            title(['total cond (g/kg), day ',num2str(tp(t)),', through y=',num2str(yplot)])
        elseif n==3; contourf(x,zt,squeeze(w(:,yplot,:))',32,'linestyle','none'); set(gca,'clim',[-5 5])
            title(['w (m/s), day ',num2str(tp(t)),', through y=',num2str(yplot)])
        elseif n==4; contourf(x,zt,squeeze(u(:,yplot,:))',32,'linestyle','none'); set(gca,'clim',[-5 5])
            title(['u (m/s), day ',num2str(tp(t)),', through y=',num2str(yplot)])
        elseif n==5; contourf(x,zt,squeeze(rv(:,yplot,:))',32,'linestyle','none'); set(gca,'clim',[0 19])
            title(['r_v (g/kg), day ',num2str(tp(t)),', through y=',num2str(yplot)])
        elseif n==6; contourf(x,zt,squeeze(RH(:,yplot,:))',32,'linestyle','none'); set(gca,'clim',[0 100])
            title(['RH (%), day ',num2str(tp(t)),', through y=',num2str(yplot)])
        end
        colorbar; ylim([0 18]); hold on; %xlim([1000 1400]);
        contour(x,zt,squeeze(totcond(:,yplot,:))','k-','levellist',[.1])
        xlabel('x (km)'); ylabel('height (km)');
    end
    
    % profiles
    figure; plot(squeeze(totcond(xplot,yplot,:)),zt,'-o'); grid on; ylabel('height (km)')
    ylim([0 15]); title(['Total cond (g/kg) through x,y=(',num2str(xplot),',',num2str(yplot),')'])
    figure; plot(squeeze(N2(xplot,yplot,:)),zt,'o-'); grid on; ylabel('height (km)')
    ylim([0 15]); title(['N^2 (s^{-1}) through x,y=(',num2str(xplot),',',num2str(yplot),')'])
    figure;plot(squeeze(theta(xplot,yplot,:)),zt,'o-'); grid on; ylabel('height (km)')
    ylim([0 15]); title(['\theta (K) through x,y=(',num2str(xplot),',',num2str(yplot),')'])
    figure; plot(squeeze(rv(xplot,yplot,:)),zt,'o-'); grid on; ylabel('height (km)')
    ylim([0 15]); title(['r_v (g/kg) through x,y=(',num2str(xplot),',',num2str(yplot),')'])
    figure; plot(squeeze(w(xplot,yplot,:)),zm,'o-'); grid on; ylabel('height (km)')
    ylim([0 15]); title(['w (m/s) through x,y=(',num2str(xplot),',',num2str(yplot),')'])
    figure; plot(squeeze(RH(xplot,yplot,:)),zt,'-o'); grid on; ylabel('height (km)')
    ylim([0 15]); title(['RH (%) through x,y=(',num2str(xplot),',',num2str(yplot),')'])
 
    clear theta rv N2 totcond w u RH
end


%% plan views read section
% max w and pcprate
if iread;
    n=0;
    for t=1:tskip:nt; n=n+1;
        w = h5read( [RAMSdir,files(t).name], '/WP', [x1 y1 z1], [nxp nyp nzp] ); % m/s
        maxw = max(w,[],3);
        maxwxt(:,n) = max(maxw,[],2);
        %wsfc(:,:,n) = w(:,:,2);
        pcprate = h5read( [RAMSdir,files(t).name], '/PCPRR', [x1 y1], [nxp nyp] )*3600.; % mm/hr
        maxpcpratext(:,n) = max(pcprate,[],2);
        totcond = (  h5read( [RAMSdir,files(t).name],'/RTP', [x1 y1 z1], [nxp nyp nzp] )...
            - h5read( [RAMSdir,files(t).name],'/RV',  [x1 y1 z1], [nxp nyp nzp] ) )*1000.; % g/kg
        maxtotcond = max(totcond,[],3);
        maxtotcondxt(:,n) = max(maxtotcond,[],2);
        thetasfc = h5read( [RAMSdir,files(t).name], '/THETA', [x1 y1 2], [nxp nyp 1] ); % K
        meanthetasfcxt(:,n) = mean(thetasfc,2);
        %rvsfc(:,:,n) = h5read( [RAMSdir,files(t).name], '/RV', [x1 y1 2], [nxp nyp 1] )*1000.; % g/kg
        %usfc(:,:,n) = h5read( [RAMSdir,files(t).name], '/UP', [x1 y1 2], [nxp nyp 1] ); % m/s
        %LWupTop = h5read( [RAMSdir,files(t).name], '/RLONTOP', [x1 y1], [nxp nyp] ); % W/m2
        LWupTop = h5read( [RAMSdir,files(t).name], '/LWUP', [x1 y1 nz-1], [nxp nyp 1] ); % W/m2
        meanOLRxt(:,n) = mean(LWupTop,2);
        LHF = h5read( [RAMSdir,files(t).name], '/SFLUX_R', [x1 y1], [nxp nyp] )*2.5e6; % W/m2
        meanLHFxt(:,n) = mean(LHF,2);
        SHF = h5read( [RAMSdir,files(t).name], '/SFLUX_T', [x1 y1], [nxp nyp] )*1004.; % W/m2
        meanSHFxt(:,n) = mean(SHF,2);
        % Precipitable water: PW = int(sfc->top){r_v*rho0*dz} [kg/m2=mm]
        rv = h5read( [RAMSdir,files(t).name], '/RV', [x1 y1 2], [nxp nyp nz-2] ); % kg/kg from z=2:nz-1
        if t==1
            % get base state density from level 2 to top-1
            cp=1004; Rd=287; fn0 = [RAMSdir,files(1).name]; % file name
            exner0 = h5read( fn0,'/PI', [x1 y1 2], [nxp nyp nz-2] )/cp; % read from 2:nz-1
            press0 = 100000.*exner0.^(cp/Rd); % pressure (Pa)
            Tvirt0 = ( exner0 .* h5read( fn0,'/THETA',[x1 y1 2],[nxp nyp nz-2] ) ) ...
                .* ( 1.+0.61*h5read( fn0,'/RV',[x1 y1 2],[nxp nyp nz-2] ) );
            rho0 = press0/Rd./Tvirt0; clear press0 Tvirt0 exner0
            % get 3D height grid
            [~,~,zzm]=meshgrid(yp,xp,zm(1:nz-1)); % zm=1:nz-1 gives PW for zt=2:nz-1
        end
        % integrate: [zm2-zm1] * rv(zt2) * rho0(zt2) + ...
        PW(:,:,n) = sum( rv .* rho0 .* diff(zzm,1,3)*1000, 3 );
        meanPWxt(:,n) = mean(PW(:,:,n),2);
        % do same for vertically integrated condensate - exclude last level in totcond
        vertintcond(:,:,n) = sum( totcond(:,:,1:end-1)/1000. .* rho0 .* diff(zzm,1,3)*1000, 3 );
        meanvintcondxt(:,n) = mean(vertintcond(:,:,n),2);
        
        if mod(tp(t),5)==0; disp(['day ',num2str(tp(t))]); end
    end
    clear w totcond rv rho0 zzm
    save([SaveDir,'planviewsXT'],'maxwxt','maxpcpratext','maxtotcondxt','meanthetasfcxt',...
        'meanOLRxt','meanLHFxt','meanSHFxt','PW','meanPWxt','meanvintcondxt')
else
    load([SaveDir,'planviewsXT']);
end


%% animation
if 0 && ny~=1
    % below line - for printing
    if exist([SaveDir,'vintcond'],'dir')~=7; mkdir([SaveDir,'vintcond']); end
    figure; set_figsize(gcf,[10 10]);%[15 5]);
    for t=1:nt%ceil(1/dt(1)):nt % once per day
        clf
        %pcolor(xp,yp,pcprate(:,:,t)'); set(gca,'clim',[0 50]); ttl='precip rate (mm/hr)';
        %pcolor(xp,yp,maxw(:,:,t)'); set(gca,'clim',[0 5]); ttl='max w (m/s)';
        %pcolor(xp,yp,maxtotcond(:,:,t)'); set(gca,'clim',[0 3]); ttl='max total cond (g/kg)';
        %pcolor(xp,yp,thetasfc(:,:,t)'); ttl='\theta at z=2 (K)'; %set(gca,'clim',[297 301]);
        %pcolor(xp,yp,rvsfc(:,:,t)'); set(gca,'clim',[13 21]); ttl='r_v at z=2 (g/kg)';
        %pcolor(xp,yp,usfc(:,:,t)'); set(gca,'clim',[-5 5]); ttl='u at z=2 (m/s)';
        %pcolor(xp,yp,wsfc(:,:,t)'); set(gca,'clim',[-.2 .2]); ttl='w at z=2 (m/s)';
        %pcolor(xp,yp,LWupTop(:,:,t)'); set(gca,'clim',[150 300]); ttl='LWup TOA (W/m^2)';
        %pcolor(xp,yp,LHF(:,:,t)'); set(gca,'clim',[0 100]); ttl='LHF (W/m^2)';
        %pcolor(xp,yp,SHF(:,:,t)'); set(gca,'clim',[0 10]); ttl='SHF (W/m^2)';
        pcolor(xp,yp,vertintcond(:,:,t)'); set(gca,'clim',[0 2]); ttl='vert int cond (mm)';
        %pcolor(xp,yp,PW(:,:,t)'); set(gca,'clim',[25 75]); ttl='PW (mm)';
        hold on; contour(xp,yp,PW(:,:,t)','w-','linewidth',1,'levellist',45)
        xlabel('x (km)'); ylabel('y (km)'); shading flat; colorbar
        title([ttl,sprintf(' - day %6.2f',tp(t))])
        %xlim([2200 2800])
        print('-dpng',[SaveDir,'vintcond/',sprintf('vintcond_day%06.2f',tp(t)),'.png'])
        pause(.25)
    end
end


%% hovmollers
for iplot=1:9
    if iplot==1; plotvar=meanPWxt; ttl='Mean (in y) PW (mm)';
    elseif iplot==2; plotvar=meanOLRxt; ttl='Mean (in y) OLR (W/m^2)';
    elseif iplot==3; plotvar=maxwxt; ttl='Max (in y and z) w (m/s)';
    elseif iplot==4; plotvar=maxpcpratext; ttl='Max (in y) precip rate (mm/hr)';
    elseif iplot==5; plotvar=maxtotcondxt; ttl='Max (in y and z) total cond (g/kg)';
    elseif iplot==6; plotvar=meanvintcondxt; ttl='Mean (in y) vert int condensate (mm)';
    elseif iplot==7; plotvar=meanthetasfcxt; ttl='Mean (in y) \theta at z=2 (K)';
    elseif iplot==8; plotvar=meanLHFxt; ttl='Mean (in y) LHF (W/m^2)';
    elseif iplot==9; plotvar=meanSHFxt; ttl='Mean (in y) SHF (W/m^2)';
    end
    figure; set_figsize(gcf,[10 15])
    contourf(xp,tp,plotvar',32,'linestyle','none'); colorbar
    xlabel('x (km)'); ylabel('time (days)'); title([ttl,' - ',ttlname2])
    if iplot==1
      disp(['Writing: ', SaveDir,'meanPW_hov']);
      print('-dpng',[SaveDir,'meanPW_hov'])
    elseif iplot==2
      disp(['Writing: ', SaveDir,'meanOLR_hov']);
      print('-dpng',[SaveDir,'meanOLR_hov'])
    end
end


%% time series

% 5-day moving averages for the smallest output interval; this means that
% if output interval changes, the 5-day moving average will actually be
% more than a 5-day moving average where dt is coarser
daymean = 5;
avgdt = ceil( daymean/min(dt) );
avgdt = avgdt + mod(avgdt+1,2); % this makes sure avgdt is an odd number

figure; set_figsize(gcf,[15 8]); plot(tp,squeeze(mean(meanPWxt,1)),'o-','linewidth',2)
hold on; plot(tp,smooth121(squeeze(mean(meanPWxt,1)),avgdt))
xlabel('time (days)'); title(['Domain-mean PW (mm) - ',ttlname2])
xlim(tp([1 end])); ylim([0 60]); grid on
disp(['Writing: ', SaveDir,'meanPW_ts']);
print('-dpng',[SaveDir,'meanPW_ts'])

load([SaveDir,'sfc-moisture.mat']) % load meanpcp variable from accumulated precip field
figure; set_figsize(gcf,[15 8]); plot(tp,meanpcp,'o-','linewidth',2)
hold on; plot(tp,smooth121(meanpcp,avgdt))
xlabel('time (days)'); title(['Domain-mean precip (mm/day) - ',ttlname2])
xlim(tp([1 end])); ylim([0 10]); grid on
disp(['Writing: ', SaveDir,'meanpcp_ts']);
print('-dpng',[SaveDir,'meanpcp_ts'])

figure; set_figsize(gcf,[15 8]); plot(tp,squeeze(max(maxwxt,[],1)),'o-','linewidth',2)
hold on; plot(tp,smooth121(squeeze(max(maxwxt,[],1)),avgdt))
xlabel('time (days)'); title(['Max w in domain (m/s) - ',ttlname2])
xlim(tp([1 end])); ylim([0 45]); grid on
disp(['Writing: ', SaveDir,'maxw_ts']);
print('-dpng',[SaveDir,'maxw_ts'])


%% cloud fraction profiles

if iread;
    n=0;
    for t=1:tskip:nt; n=n+1;
        totcond = (  h5read( [RAMSdir,files(t).name], '/RTP', [x1 y1 z1], [nxp nyp nzp] ) ...
                   - h5read( [RAMSdir,files(t).name], '/RV',  [x1 y1 z1], [nxp nyp nzp] ) ) *1000.; % g/kg
        CF = zeros(size(totcond)); CF(totcond>=0.01)=1; % cloudy grid points: where total cond >= 0.01 g/kg
        CFprof(:,n) = sum(sum(CF,1),2)/(nxp*nyp)*100; % cloud fraction in %
        CFxz(:,:,n) = sum(CF,2)/nyp*100; % cloud fraction cross section
        meanTCxz(:,:,n) = mean(totcond,2);
        if mod(tp(t),5)==0; disp(['day ',num2str(tp(t))]); end
    end
    save([SaveDir,'CF'],'CFprof','CFxz','meanTCxz')
    clear CF totcond
else
    load([SaveDir,'CF']);
end

figure; set_figsize(gcf,[15 8])
contourf(tp,zp,CFprof,32,'linestyle','none'); colorbar; set(gca,'clim',[0 50])
title(['Cloud fraction (%) - ',ttlname2])
xlabel('time (days)'); ylabel('height (km)')
disp(['Writing: ', SaveDir,'CF_hov']);
print('-dpng',[SaveDir,'CF_hov'])

figure; set_figsize(gcf,[15 8])
contourf(tp,zp,squeeze(mean(meanTCxz,1)),32,'linestyle','none'); colorbar; %set(gca,'clim',[0 5])
title(['Mean total cond (g/kg) - ',ttlname2])
xlabel('time (days)'); ylabel('height (km)')
disp(['Writing: ', SaveDir,'meanTC_hov']);
print('-dpng',[SaveDir,'meanTC_hov'])

% mean over last 10 days
if tp(end)-10>0 % there are enough days
    itm10days=find(tp(1:tskip:nt)>=tp(end)-10,1,'first');
else itm10days=1; % index 1
end
figure; set_figsize(gcf,[10 10])
plot(mean(CFprof(:,itm10days:end),2),zp); title(['Mean cloud fraction last 10 days (%) - ',ttlname2])
set(gca,'ylim',zp([1 end]))
ylabel('height (km)')
grid on
disp(['Writing: ', SaveDir,'meanCF_prof']);
print('-dpng',[SaveDir,'meanCF_prof'])

% mean over last 10 days
figure; set_figsize(gcf,[10 10])
plot(squeeze(mean(mean(meanTCxz(:,:,itm10days:end),1),3)),zp); title(['Mean total cond last 10 days (g/kg) - ',ttlname2])
set(gca,'ylim',zp([1 end]))
ylabel('height (km)')
grid on
disp(['Writing: ', SaveDir,'meanTC_prof']);
print('-dpng',[SaveDir,'meanTC_prof'])

%figure; set_figsize(gcf,[15 8]); n=0;
%for t=1:tskip:nt; n=n+1;
%    clf
%    pcolor(xp,zp,CFxz(:,:,n)'); set(gca,'clim',[0 70]); title(['Cloud fraction (%) - day ',num2str(tp(t))])
%    xlabel('x (km)'); ylabel('z (km)'); shading flat; colorbar
%    pause(.5)
%end


%% hydrometeors

if 0
mvars = {'RCP';'RRP';'RPP';'RSP';'RAP';'RGP';'RHP';'RV'};
mvarnames = {'cloud';'rain';'pristine ice';'snow';'aggregates';'graupel';'hail';'vapor'};

if iread;
    n=0;
    mvarxz=zeros(nxp,nzp,length(1:tskip:nt),length(mvars));
    for v = 1:length(mvars); n=0;
        for t=1:tskip:nt; n=n+1;
            mvar = h5read( [RAMSdir,files(t).name], ['/',mvars{v}], [x1 y1 z1], [nxp nyp nzp] );
            mvarxz(:,:,n,v) = mean(mvar,2)*1000.; % cross section (g/kg)
        end
        disp(mvars{v})
    end
    save([SaveDir,'MicroMRs'],'mvarxz')
else
    load([SaveDir,'MicroMRs']);
end

for v=1:length(mvars)
    figure; set_figsize(gcf,[15 8])
    contourf(tp(1:tskip:nt),zp,squeeze(mean(mvarxz(:,:,:,v),1)),32,'linestyle','none'); colorbar; %set(gca,'clim',[0 5])
    title(['Mean ',mvarnames{v},' (g/kg) - ',ttlname2])
    xlabel('time (days)'); ylabel('height (km)')
    print('-dpng',[SaveDir,'mean',mvars{v},'_hov'])
end

% mean over last 10 days
if tp(end)-10>0 % there are enough days
    itm10days=find(tp(1:tskip:nt)>=tp(end)-10,1,'first');
else itm10days=1; % index 1
end
figure; set_figsize(gcf,[10 10]); %load_cmaps; 
set(gcf,'DefaultAxesColorOrder',parula(round(linspace(1,64,length(mvars))),:)); hold all
plot(squeeze(mean(mean(sum(mvarxz(:,:,itm10days:end,1:end-1),4),1),3)),zp); % total cond
for v=1:length(mvars)-1 % exclude water vapor
    plot(squeeze(mean(mean(mvarxz(:,:,itm10days:end,v),1),3)),zp); 
end
legend({'total cond',mvarnames{1:end-1}})
title(['Mean mixing ratios last 10 days (g/kg) - ',ttlname2])
set(gca,'ylim',zp([1 end]))
ylabel('height (km)')
grid on; box on
print('-dpng',[SaveDir,'meanMRs_prof'])

end

%% mean cross sections in equilibrium
if nx>5*ny % only do this section for channel domain...
    it1=find(tp>=45,1); it2=nt;
    if isempty(it1); it1=nt+1; end
    if iread
        n=0;
        for t=1:nt; n=n+1;
            tempK = h5read( [RAMSdir,files(t).name],'/THETA', [x1 y1 z1], [nxp nyp nzp] ); % potential temp K
            rv = h5read( [RAMSdir,files(t).name],'/RV', [x1 y1 z1], [nxp nyp nzp] )*1000.; % g/kg
            rvxc(:,:,n) = mean(rv,2);
            press = h5read( [RAMSdir,files(t).name],'/PI', [x1 y1 z1], [nxp nyp nzp] )/1004; % pi dimensionless
            tempK = tempK.*press; % (K) [here press is still pi]
            press = 1000.*press.^(1004/287); % mb
            RH = rsat( tempK, press, rv, 'var','rhw' ); % pct
            RHxc(:,:,n) = mean(RH,2);
            
            totcond = ( h5read( [RAMSdir,files(t).name],'/RTP', [x1 y1 z1], [nxp nyp nzp] )...
                -h5read( [RAMSdir,files(t).name],'/RV', [x1 y1 z1], [nxp nyp nzp] ) )*1000.;
            totcondxc(:,:,n) = mean(totcond,2);
            
            w = h5read( [RAMSdir,files(t).name],'/WP', [x1 y1 z1], [nxp nyp nzp] );
            wxc(:,:,n) = mean(w,2);
            % u: use variable w
            w = h5read( [RAMSdir,files(t).name],'/UP', [x1 y1 z1], [nxp nyp nzp] );
            uxc(:,:,n) = mean(w,2);
            
            if mod(tp(t),5)==0; disp(['day ',num2str(tp(t))]); end
        end
        save([SaveDir,'MeanxcEquil'],'rvxc','RHxc','totcondxc','wxc','uxc')
        clear tempK rv press RH totcond w
    else
        load([SaveDir,'MeanxcEquil'])
    end
    
    if exist([SaveDir,'equilXCs'],'dir')~=7; mkdir([SaveDir,'equilXCs']); end
    %load_cmaps;
    
    for ivar=1:5
        if ivar==1; pvar=rvxc; ttl='r_v (g/kg)'; pttl='rv'; pclim=[0 20];
        elseif ivar==2; pvar=RHxc; ttl='RH (%)'; pttl='RH'; pclim=[0 100];
        elseif ivar==3; pvar=totcondxc; ttl='total cond (g/kg)'; pttl='totcond'; pclim=[0 1];
        elseif ivar==4; pvar=wxc; ttl='w (m/s)'; pttl='w'; pclim=[-.1 .1];
        elseif ivar==5; pvar=uxc; ttl='u (m/s)'; pttl='u'; pclim=[-5 5];
        end
        if exist([SaveDir,'equilXCs/',pttl],'dir')~=7; mkdir([SaveDir,'equilXCs/',pttl]); end
        figure; set_figsize(gcf,[15 5]);%[10 10]);%
        n=0;
        for t=1:nt+(it2>it1); n=n+1;
            clf
            if t<=nt; pcolor(xp,zp,pvar(:,:,n)')
                title(['mean (in y) ',ttl,sprintf(' - day %6.2f',tp(t))])
            else % this will only happen if run is far enough
                pcolor(xp,zp,mean(pvar(:,:,it1:it2),3)'); title(['mean ',ttl])
            end
            shading flat; colorbar;
            if t==nt+1
                if ivar==3; set(gca,'clim',pclim/5);
                elseif ivar==4; set(gca,'clim',pclim/2);
                end
            else set(gca,'clim',pclim);
            end
            if any(ivar==[4 5]); set(gcf,'colormap',redblu); ylim([0 20]);
            else ylim([0 15]);
            end
            xlabel('x (km)'); ylabel('height (km)');
            if t<=nt; print('-dpng',[SaveDir,'equilXCs/',pttl,'/mean',pttl,sprintf('_day%06.2f',tp(t)),'.png'])
            else print('-dpng',[SaveDir,'equilXCs/mean',pttl,sprintf('_day%06.2f-%06.2f',tp(it1),tp(it2)),'.png'])
            end
            pause(.5)
        end
    end
end

disp('Done:');

%% Copy .png files to google drive

if SaveGoogleDrive==1

    if ~exist(GoogleDrivePath,'dir'); mkdir(GoogleDrivePath); end
    
    copyfile([SaveDir,'*.png'],GoogleDrivePath)
    copyfile([SaveDir,'vintcond'],[GoogleDrivePath,'vintcond'])
    
    % copy some specific files and also rename them to include the run name
    copyfile([SaveDir,'CF_hov.png'],[GoogleDrivePathBase,'BasicPlots_Compare/CF_hov_',ttlname,'.png'])
    copyfile([SaveDir,'meanMRs_prof.png'],[GoogleDrivePathBase,'BasicPlots_Compare/meanMRs_prof_',ttlname,'.png'])
    copyfile([SaveDir,'avgrv_prof.png'],[GoogleDrivePathBase,'BasicPlots_Compare/avgrv_prof_',ttlname,'.png'])
    copyfile([SaveDir,'avgtemp_prof.png'],[GoogleDrivePathBase,'BasicPlots_Compare/avgtemp_prof_',ttlname,'.png'])
    copyfile([SaveDir,'avgRH_prof.png'],[GoogleDrivePathBase,'BasicPlots_Compare/avgRH_prof_',ttlname,'.png'])
    
end
